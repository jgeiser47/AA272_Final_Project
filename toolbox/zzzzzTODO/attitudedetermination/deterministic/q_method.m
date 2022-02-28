%==========================================================================
%
% q_method  Statistical attitude determination.
%
%   q_body = q_method(M,U,w)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   M       (3xN) measurement matrix
%   U       (3xN) reference matrix
%   w       (Nx1) weight vector
%
% --------
% OUTPUTS:
% --------
%   q_body	(4x1) body attitude quaternion (parameterizes attitude of body
%                 frame w.r.t. ECI frame)
%
% -----
% NOTE:
% -----
%   --> Each column of "M" is a unit vector measurement resolved in the
%       measurement (body) frame
%   --> Each column of "U" is the corresponding ground truth reference unit
%       vector resolved in the reference (ECI) frame.
%   --> Each element of "w" stores the relative weight (i.e. importance) of
%       the corresponding measurement.
%
%==========================================================================
function q_body = q_method(M,U,w)
    
    % determines number of measurements
    N = size(M,2);
    
    % preallocates W and V matrices
    W = zeros(3,N);
    V = zeros(3,N);
    
    % assembles W and V matrices
    for i = 1:N
        W(:,i) = sqrt(w(i))*M(:,i);
        V(:,i) = sqrt(w(i))*U(:,i);
    end
    
    % assembles matrices and calculates necessary quantities
    B = W*V';
    S = B'+B;
    Z = [B(2,3)-B(3,2);B(3,1)-B(1,3);B(1,2)-B(2,1)];
    sigma = trace(B);
    K = [S-eye(3)*sigma,Z;Z',sigma];
    
    % eigenvector and eigenvalue matrices of K
    [Q,Lambda] = eig(K);

    % deals with error where Simulink assumes eigenvalues can be complex
    Q = real(Q);
    Lambda = real(Lambda);
    
    % vector of eigenvalues
    lambda = diag(Lambda);
    
    % index of maximum eigenvalue
    [~,imax] = max(lambda);
    
    % extracts the attitude quaternion from the eigenvector matrix
    q_body = Q(:,imax);
    
end