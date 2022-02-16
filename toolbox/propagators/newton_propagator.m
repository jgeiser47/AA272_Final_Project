%==========================================================================
%
% newton_propagator  Time derivative of a satellite's ECI state vector.
%
%   [dXdt,extra] = newton_propagator(t,x,prop,sat)
%
% Author: Tamas Kis
% Last Update: 2022-02-07
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   t           - (1×1 double) simulation time [s]
%   x           - (6×1 double) orbital state vector
%                   --> (1-3): r_eci - position resolved in ECI frame [m]
%                   --> (4-6): v_eci - inertial velocity resolved in ECI
%                                      frame [m/s]
%   prop        - (struct) propagator settings
%   sat  - (struct) satellite parameters
%               
%
% -------
% OUTPUT:
% -------
%   dx	        - (6×1 double) orbital state vector derivative
%                   --> (1-3): v_eci - inertial velocity resolved in ECI
%                                      frame [m/s]
%                   --> (4-6): a_eci - inertial acceleration resolved in
%                                      ECI frame [m/s^2]
%   extra       - (struct) extra parameters to return
%       • a_eci    - (3×1 double) inertial acceleration resolved in ECI 
%                    frame [m/s^2]
%       • r_ecef   - (3×1 double) position resolved in ECEF frame [m]
%       • v_ecef   - (3×1 double) ECEF velocity resolved in ECEF frame 
%                    [m/s]
%       • lat      - (1×1 double) geodetic latitude [deg]
%       • lon      - (1×1 double) geodetic longitude [deg]
%       • h        - (1×1 double) geodetic altitude [deg]
%       • rho      - (1×1 double) atmospheric density [kg/m^3]
%
%==========================================================================
function [dXdt,extra] = newton_propagator(t,X,prop,sat)

    % TODO: update DUT1/DAT

    % extracts ECI position [m] and inertial velocity [m/s] 
    r_eci = X(1:3);
    v_eci = X(4:6);

    % -----------------------------
    % Timing and Earth orientation.
    % -----------------------------

    % time scales [MJD]
    [MJD_GPS,MJD_TAI,MJD_TT,MJD_UT1,MJD_UTC] = time_scales(t,prop.t0);

    % Gregorian date of UTC [y,mo,d,h,m,s]
    cal_UTC = mjd2cal(MJD_UTC);

    % Earth orientation parameters for IAU2006/2000 CIO based theory
    [dX,dY,xp,yp,LOD] = eop_iau06(MJD_UT1,prop.data.eop);

    % rotation matrix (GCRF --> ITRF) and Earth angular velocity resolved
    % in the ITRF [rad/s] from IAU2006/2000 CIO based theory
    [R_ecef2eci,R_eci2ecef,w_eci] = iau06(MJD_UT1,MJD_TT,xp,yp,dX,dY,...
        LOD,prop.data.XYs_iau06);
    
    % ----------------------------
    % Other state representations.
    % ----------------------------

    % ECEF position [m] and velocity [m/s]
    [r_ecef,v_ecef] = eci2ecef(r_eci,v_eci,w_eci,R_eci2ecef);
    
    % geodetic latitude [rad], longitude [rad], and altitude [m]
    [lat,lon,h] = ecef2geod(r_ecef);

    % ------------
    % Environment.
    % ------------

    % Sun position resolved in ECI frame [m]
    r_sun_eci = sun_position(MJD_UT1,MJD_TT);

    % determines if satellite is in eclipse
    in_eclipse = eclipse(r_eci,r_sun_eci);

    % atmospheric density [kg/m^3]
    if strcmpi(prop.models.density,'Harris-Priester')
        rho = density_harris_priester(r_eci,r_ecef,r_sun_eci);
    elseif strcmpi(prop.models.density,'Exponential')
        rho = density_exponential(r_ecef);
    elseif strcmpi(prop.models.density,'Jacchia-Roberts')
        rho = density_jacchia_roberts(r_eci,r_ecef,r_sun_eci,...
            R_eci2ecef,MJD_TT);
    end

    % --------------
    % Accelerations.
    % --------------

    % acceleration due to gravity resolved in ECI frame [m/s^2]
    g_eci = R_ecef2eci*gravity(r_ecef,prop.models.grav_mu,...
        prop.models.grav_R,prop.data.C,prop.data.S,prop.models.grav_N);
    %g_eci = -(MU_EARTH/inorm(r_eci)^3)*r_eci+R_ecef2eci*J2_ecef(r_ecef);

    % perturbing acceleration due to atmospheric drag resolved in ECI frame
    % [m/s^2]
    if prop.perturb.drag
        fD_eci = drag(v_ecef,R_ecef2eci,rho,sat.B);
    else
        fD_eci = [0;0;0];
    end
    
    % perturbing acceleration due to solar radiation pressure resolved
    % in ECI frame[m/s^2]
    if prop.perturb.srp && ~in_eclipse
        f_srp_eci = srp(r_eci,r_sun_eci,sat.CR,sat.Asrp,sat.M);
    else
        f_srp_eci = [0;0;0];
    end
    
    % perturbing acceleration due to general relativity resolved in ECI
    % frame
    if prop.perturb.relativity
        f_rel_eci = relativity(r_eci,v_eci);
    else
        f_rel_eci = [0;0;0];
    end

    % total inertial acceleration resolved in ECI frame [m/s^2]
    a_eci = g_eci+fD_eci+f_srp_eci+f_rel_eci;
    
    % ------------------------
    % State vector derivative.
    % ------------------------
    
    dXdt = [v_eci;a_eci];

    % -------------------------
    % Storing extra parameters.
    % -------------------------

    extra.a_eci = a_eci;
    extra.r_ecef = r_ecef;
    extra.v_ecef = v_ecef;
    extra.lat = lat;
    extra.lon = lon;
    extra.h = h;

    extra.rho = rho;
    extra.in_eclipse = in_eclipse;
    extra.r_sun_eci = r_sun_eci;

    extra.fD_eci = fD_eci;
    extra.f_srp_eci = f_srp_eci;
    extra.f_rel_eci = f_rel_eci;
    extra.g_eci = g_eci;

    extra.MJD_GPS = MJD_GPS;
    extra.MJD_TAI = MJD_TAI;
    extra.MJD_TT = MJD_TT;
    extra.MJD_UT1 = MJD_UT1;
    extra.MJD_UTC = MJD_UTC;
    extra.cal_UTC = cal_UTC;
    
end