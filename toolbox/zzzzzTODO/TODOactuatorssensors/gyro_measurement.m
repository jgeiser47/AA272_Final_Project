%==========================================================================
%
% gyro_measurement  Simulates gyroscope measurements.
%
%   wm_body = gyro_measurement(w,K_gyro,b_gyro,Sigma_gyro,R_body2prin)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   w           (3x1) [rad/s] true angular velocity vector resolved in 
%                             principal frame 
%   K_gyro      (3x3) scale factor matrix
%   b_gyro      (3x1) [rad/s] drift (bias)
%   Sigma_gyro  (3x3) [rad^2/s^2] noise covariance
%   R_body2prin (3x3) rotation matrix (body --> principal)
%
% --------
% OUTPUTS:
% --------
%   wm_body     (3x1) [rad/s] angular velocity measurement resolved in body
%                             frame
%
%==========================================================================
function wm_body = gyro_measurement(w,K_gyro,b_gyro,Sigma_gyro,R_body2prin)
    
    % resolves true angular velocity in body frame [rad/s]
    w_body = R_body2prin'*w;
    
    % sets up gaussian_random_sample MATLAB function for use in Simulink
    coder.extrinsic('gaussian_random_sample');
    
    % initializes noise to avoid Simulink errors
    n = zeros(3,1);
    
    % samples noise [rad/s]
    n = gaussian_random_sample(zeros(3,1),Sigma_gyro);
    
    % simulates measured body frame angular velocity [rad/s]
    wm_body = K_gyro*w_body+b_gyro+n;
    
end