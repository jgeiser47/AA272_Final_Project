%==========================================================================
%
% dynamics_model  Dynamics model TODO. 
%
%   dxdt = dynamics_model(t,x,filter,sat)
%
% Author: Tamas Kis
% Last Update: 2022-02-21
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   t       - (1×1 double) simulation time [s]
%   x       - (6×1 double) state vector
%               --> 1-3. r_eci - position resolved in ECI frame [m]
%               --> 4-6. v_eci - inertial velocity resolved in ECI frame 
%                                [m/s]
%   filter  - (1×1 struct) data for orbit propagation/determination (see 
%             "initialize_filter" for full definition
%   sat     - (1×1 struct) satellite parameters TODO   
%
% -------
% OUTPUT:
% -------
%   dxdt    - (6×1 double) state vector derivative
%               --> 1-3. v_eci - inertial velocity resolved in ECI frame
%                                [m/s]
%               --> 4-6. a_eci - inertial acceleration resolved in ECI
%                                frame [m/s^2]
%
%==========================================================================
function dxdt = dynamics_model(t,x,filter,sat)
    
    % extracts ECI position [m] and inertial velocity [m/s] 
    r_eci = x(1:3);
    v_eci = x(4:6);
    
    % -----------------------------
    % Timing and Earth orientation.
    % -----------------------------

    % time scales [MJD]
    [~,~,MJD_TT,MJD_UT1] = time_scales(t,filter.t0);

    % Earth orientation parameters for IAU2006/2000 CIO based theory
    [xp,yp,dX,dY,LOD] = eop_iau06(MJD_UT1,filter.data.eop);

    % rotation matrix (GCRF --> ITRF) and Earth angular velocity resolved
    % in the ITRF [rad/s] from IAU2006/2000 CIO based theory
    [R_ecef2eci,R_eci2ecef,w_eci] = iau06(MJD_UT1,MJD_TT,xp,yp,dX,dY,...
        LOD,filter.data.XYs_iau06);
    
    % ----------------------------
    % Other state representations.
    % ----------------------------

    % ECEF position [m]
    r_ecef = eci2ecef(r_eci,v_eci,w_eci,R_eci2ecef);

    % ------------------------------------------
    % Accelerations (all resolved in ECI frame).
    % ------------------------------------------
    
    % acceleration due to gravity [m/s^2]
    g_eci = R_ecef2eci*gravity(r_ecef,filter.data.C,filter.data.S,2);

    % total inertial acceleration resolved in ECI frame [m/s^2]
    a_eci = g_eci;
    
    % ------------------------
    % State vector derivative.
    % ------------------------
    
    dxdt = [v_eci;a_eci];

end