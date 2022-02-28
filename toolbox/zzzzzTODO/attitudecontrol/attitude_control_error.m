%==========================================================================
%
% attitude_control_error  Determines the attitude control error
% parameterized as 3-2-1 Euler angles.
%
%   [dphi_ctrl,dtheta_ctrl,dpsi_ctrl] = attitude_control_error(A,A_nom)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% INPUTS:
%   A       actual (can be estimated or ground-truth) attitude matrix
%           (3 x 3 x n)
%   A_nom   nominal attitude matrix (3 x 3 x n)
%
% OUTPUTS:
%   dphi_ctrl       attitude ctrl error as 3-2-1 Euler angle (n x 1) [rad]
%   dtheta_ctrl     attitude ctrl error as 3-2-1 Euler angle (n x 1) [rad]
%   dpsi_ctrl       attitude ctrl error as 3-2-1 Euler angle (n x 1) [rad]
%
%==========================================================================
function [dphi_ctrl,dtheta_ctrl,dpsi_ctrl] = attitude_control_error(...
    A,A_nom)

    % determines number of attitude matrices
    n = size(A_nom,3);
    
    % preallocates arrays
    dphi_ctrl = zeros(n,1);
    dtheta_ctrl = zeros(n,1);
    dpsi_ctrl = zeros(n,1);
    
    % attitude control error parameterized as 3-2-1 Euler angles for
    % each set of estimated and ground-truth attitudes
    for i = 1:n
        dA_ctrl = A(:,:,i)*A_nom(:,:,i)';
        [dphi_ctrl(i),dtheta_ctrl(i),dpsi_ctrl(i)] = dcm2eul321(dA_ctrl);
    end
    
    % unwraps to find "cumulative" error
    dphi_ctrl = unwrap(dphi_ctrl);
 	dtheta_ctrl = unwrap(dtheta_ctrl);
 	dpsi_ctrl = unwrap(dpsi_ctrl);
    
end