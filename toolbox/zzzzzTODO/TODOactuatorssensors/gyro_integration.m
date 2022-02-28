%==========================================================================
%
% gyro_integration  Integrates the body frame angular velocity measured by
% a gyroscope to obtain the attitude quaternion.
%
%   q_new = gyro_integration(q_old,wm_body,Tg)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   q_old   (4x1) body attitude quaternion (parameterizes attitude of body
%                 frame w.r.t. ECI frame) at current sampling time
%   wm_body (3x1) [rad/s] angular velocity measurement resolved in body
%                         frame
%   T       (1x1) [s] gyroscope sampling period
%
% --------
% OUTPUTS:
% --------
%   q_body	(4x1) body attitude quaternion at next sampling time 
%
%==========================================================================
function q_new = gyro_integration(q_old,wm_body,Tg)
    
    % unpacks angular velocity vector [rad/s]
    wxb = wm_body(1);
    wyb = wm_body(2);
    wzb = wm_body(3);
    
    % angular velocity matrix for quaternion propagation [rad/s]
    Om = [ 0      wzb   -wyb    wxb; 
          -wzb    0      wxb    wyb;
           wyb   -wxb    0      wzb;
          -wxb   -wyb   -wzb    0];
    
    % angular velocity magnitude [rad/s]
    w = norm(wm_body);
    
    % propagates quaternion
    q_new = (cos(w*Tg/2)*eye(4)+(1/w)*sin(w*Tg/2)*Om)*q_old;
    
    % renormalizes quaternion
    q_new = q_new/norm(q_new);

end