%==========================================================================
%
% attitude_estimation_error  Determines the attitude estimation error.
%
%   [dphi_est,dtheta_est,dpsi_est] = attitude_estimation_error(A_est,A)
%   [dq1_est,dq2_est,dq3_est,dq4_est] = attitude_estimation_error(A_est,A)
%
% Author: Tamas Kis
% Last Update: 2022-01-05
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   A_est       - (3×3×n double) estimated attitude matrix
%   A           - (3×3×n double) ground-truth attitude matrix
%   param       - (OPTIONAL) (char) parameterization --> 'eul321' for 3-2-1
%                 Euler angles or 'quat' for quaternions
%
% ----------------------
% OUTPUT (for 'eul321'):
% ----------------------
%   dphi_est    - (n×1 double) 1st 3-2-1 Euler angle of attitude estimation
%                 error [rad]
%   dtheta_est  - (n×1 double) 2nd 3-2-1 Euler angle of attitude estimation
%                 error [rad]
%   dpsi_est    - (n×1 double) 3rd 3-2-1 Euler angle of attitude estimation
%                 error [rad]
%
% --------------------
% OUTPUT (for 'quat'):
% --------------------
%   dq1_est     - (n×1 double) 1st quaternion element of attitude 
%                 estimation error
%   dq2_est     - (n×1 double) 2nd quaternion element of attitude 
%                 estimation error
%   dq3_est     - (n×1 double) 3rd quaternion element of attitude 
%                 estimation error
%   dq4_est     - (n×1 double) 4th quaternion element of attitude 
%                 estimation error
%
% -----
% NOTE:
% -----
%   --> n - number of attitude matrices input
%
%==========================================================================
function varargout = attitude_estimation_error(A_est,A,parameterization)

    % sets default parameterization to 3-2-1 Euler angles
    if (nargin < 3) || isempty(parameterization)
        parameterization = 'eul321';
    end
    
    % determines number of attitude matrices
    n = size(A,3);
    
    % -------------------------------------------------------------
    % Attitude estimation error using Euler 3-2-1 parameterization.
    % -------------------------------------------------------------

    if strcmp(parameterization,'eul321')
        
        % preallocates arrays
        dphi_est = zeros(n,1);
        dtheta_est = zeros(n,1);
        dpsi_est = zeros(n,1);

        % attitude estimation error parameterized as 3-2-1 Euler angles for
        % each set of estimated and ground-truth attitudes
        for i = 1:n
            dA_est = A_est(:,:,i)*A(:,:,i)';
            [dphi_est(i),dtheta_est(i),dpsi_est(i)] = dcm2eul321(dA_est);
        end

        % unwraps to find "cumulative" error
        dphi_est = unwrap(dphi_est);
        dtheta_est = unwrap(dtheta_est);
        dpsi_est = unwrap(dpsi_est);
        
        % outputs results
        varargout = {dphi_est,dtheta_est,dpsi_est};
        
    % ------------------------------------------------------------
    % Attitude estimation error using quaternion parameterization.
    % ------------------------------------------------------------

    elseif strcmp(parameterization,'quat')
        
        % preallocates arrays
        dq1_est = zeros(n,1);
        dq2_est = zeros(n,1);
        dq3_est = zeros(n,1);
        dq4_est = zeros(n,1);

        % attitude estimation error parameterized as quaternions for each
        % set of estimated and ground-truth attitudes
        for i = 1:n
            dA_est = A_est(:,:,i)*A(:,:,i)';
            dq_est = dcm2quat(dA_est);
            dq1_est(i) = dq_est(1);
            dq2_est(i) = dq_est(1);
            dq3_est(i) = dq_est(1);
            dq4_est(i) = dq_est(1);
        end
        
        % outputs results
        varargout = {dq1_est,dq2_est,dq3_est,dq4_est};
        
    end
    
end