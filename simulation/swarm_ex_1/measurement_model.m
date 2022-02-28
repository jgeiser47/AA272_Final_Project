%==========================================================================
%
% measurement_model  Measurement model TODO. 
%
%   y = measurement_model(t,x,filter,sat)
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
%   y       - (6×1 double) measurement
%               --> 1-3. r_ecef - inertial velocity resolved in ECI frame
%                                 [m/s]
%               --> 4-6. v_ecef - ECEF velocity resolved in ECEF frame 
%                                 [m/s]
%
%==========================================================================
function y = measurement_model(t,x,filter,sat)
    
    % extracts ECI position [m] and inertial velocity [m/s] 
    r_eci = x(1:3);
    v_eci = x(4:6);

    % time scales [MJD]
    [~,~,MJD_TT,MJD_UT1] = time_scales(t,filter.t0);
    
    % Earth orientation parameters for IAU2006/2000 CIO based theory
    [xp,yp,dX,dY,LOD] = eop_iau06(MJD_UT1,filter.data.eop);
    
    % rotation matrix (GCRF --> ITRF) and Earth angular velocity resolved
    % in the ITRF [rad/s] from IAU2006/2000 CIO based theory
    [~,R_eci2ecef,w_eci] = iau06(MJD_UT1,MJD_TT,xp,yp,dX,dY,...
        LOD,filter.data.XYs_iau06);

    % position [m] and ECEF velocity [m/s] resolved in ECEF frame
    [r_ecef,v_ecef] = eci2ecef(r_eci,v_eci,w_eci,R_eci2ecef);

    % measurement
    y = [r_ecef;
         v_ecef];

end