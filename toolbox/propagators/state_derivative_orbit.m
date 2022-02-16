%==========================================================================
%
% state_derivative_orbit  Derivative of a satellite's orbital state vector.
%
%   [dx,extra] = state_derivative_orbit(t,x,prop,sat)
%
% Author: Tamas Kis
% Last Update: 2022-01-17
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
%   spacecraft  - (struct) spacecraft parameters
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
function [dx,extra] = state_derivative_orbit(t,x,prop,spacecraft)

    % TODO: add control input

    % extracts ECI position [m] and inertial velocity [m/s] 
    r_eci = x(1:3);
    v_eci = x(4:6);

    % -----------------------------
    % Timing and Earth orientation.
    % -----------------------------

    % time parameters
    [GPS_cal,GPS_MJD,GPS_ws,TT_MJD,UT1_MJD,UTC_cal,UTC_MJD] =...
        simtime2realtime(t,prop.epoch);

    % Greenwich mean sidereal time [rad]
    GMST = UT1toGMST(UT1_MJD);

    % rotation matrix from ECEF frame to ECI frame
    R_ecef2eci = ecef2eci_matrix(GMST);
    
    % ----------------------------
    % Other state representations.
    % ----------------------------

    % ECEF position [m] and velocity [m/s]
    [r_ecef,v_ecef] = eci2ecef(r_eci,v_eci,GMST);

    % geodetic latitude [rad], longitude [rad], and altitude [m]
    [lat,lon,h] = ecef2geod(r_ecef);

    % ------------
    % Environment.
    % ------------

    % atmospheric density [kg/m^3]
    rho = density_harris_priester(r_eci,r_ecef,UT1_MJD);

    % --------------
    % Accelerations.
    % --------------

    % acceleration due to atmospheric drag resolved in ECI frame [m/s^2]
    fD_eci = accel_drag(r_eci,v_eci,rho,spacecraft.B);

    % acceleration due to gravity resolved in ECI [m/s^2]
    %g_eci = accel_gravity(r_eci,r_ecef);

    % acceleration due to gravity resolved in ECEF frame [m/s^2]
    a_ecef = gravity(r_ecef,mu,R,C,S);

    % acceleration due to gravity resolved in ECI frame [m/s^2]
    a_eci = R_ecef2eci*a_ecef;

    % total inertial acceleration resolved in ECI frame [m/s^2]
    a_eci = a_eci+fD_eci;
    
    % packages state vector derivative
    dx = [v_eci;a_eci];

    % stores extra parameters
    extra.a_eci = a_eci;
    extra.r_ecef = r_ecef;
    extra.v_ecef = v_ecef;
    extra.lat = lat;
    extra.lon = lon;
    extra.h = h;
    extra.rho = rho;
    
end