%==========================================================================
%
% eci2rtn  Deputy relative state in chief's RTN frame from the ECI states
% of the chief and the deputy.
%
%   dX_rtn = eci2rtn(Xc_eci,Xd_eci)
%
% Author: Tamas Kis
% Last Update: 2022-02-21
%
% REFERENCES:
%   [1] Curtis, "Orbital Mechanics for Engineering Students", 3rd Ed.,
%       Eqs. (7.5), (7.7), (7.8), and (7.12) (pp. 369-370)
%   [2] D'Amico, "Orbit Perturbations, Gauss Variational Equations", AA 
%       279A Lecture 9, Personal Notes (p. 4)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   Xc_eci	- (6×1 double) chief ECI state
%               1-3. position resolved in ECI frame
%               4-6. inertial velocity resolved in ECI frame
%   Xd_eci  - (6×1 double) deputy ECI state
%               1-3. position resolved in ECI frame
%               4-6. inertial velocity resolved in ECI frame
%
% -------
% OUTPUT:
% -------
%   dX_rtn 	- (6×1 double) deputy relative state in chief's RTN frame
%               1-3. relative position resolved in chief RTN frame
%               4-6. relative velocity resolved in chief RTN frame
%
% -----
% NOTE:
% -----
%   --> The states can be input in any units, but they MUST be consistent. 
%       The output state will be in the same units.
%
%==========================================================================
function dx_rtn = eci2rtn(xc_eci,xd_eci)

    % extracts ECI positions and velocities from ECI states
    rc_eci = xc_eci(1:3);
    vc_eci = xc_eci(4:6);
    rd_eci = xd_eci(1:3);
    vd_eci = xd_eci(4:6);
    
    % angular velocity of RTN frame w.r.t. ECI frame, resolved in ECI frame
    % [rad/s]
    w_rtn_eci = cross(rc_eci,vc_eci)/norm(rc_eci)^2;
    
    % relative position and velocity resolved in ECI frame
    dr_eci = rd_eci-rc_eci;
    dv_eci = vd_eci-vc_eci-cross(w_rtn_eci,dr_eci);

    % R, T, and N basis vectors resolved in ECI frame
    R_eci = rc_eci/norm(rc_eci);
    N_eci = cross(rc_eci,vc_eci)/norm(cross(rc_eci,vc_eci));
    T_eci = cross(N_eci,R_eci);
    
    % rotation matrix from ECI frame to RTN frame
    R_eci2rtn = [R_eci,T_eci,N_eci]';
    
    % relative position and velocity resolved in RTN frame
    dr_rtn = R_eci2rtn*dr_eci;
    dv_rtn = R_eci2rtn*dv_eci;
    
    % assembles RTN state of deputy
    dx_rtn = [dr_rtn;dv_rtn];
    
end