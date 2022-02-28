%==========================================================================
%
% magnetometer_measurement  Simulates magnetometer measurements.
%
%   Bm_body = magnetometer_measurement(B_eci,b_mag,Sigma_mag,A,R_body2prin)
%
% Copyright (c) 2021 Luke Neise, Tamas Kis
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   B_eci       (3x1) [T] true magnetic field vector resolved in ECI frame
%   b_mag       (3x1) [T] bias
%   Sigma_mag   (3x3) [T^2] noise covariance
%   A           (3x3) attitude matrix (ECI --> principal)
%   R_body2prin (3x3) rotation matrix (body --> principal)
%
% --------
% OUTPUTS:
% --------
%   Bm_body     (3x1) [T] magnetic field vector measurement resolved in 
%               body frame
%
%==========================================================================
function Bm_body = magnetometer_measurement(B_eci,b_mag,Sigma_mag,A,...
    R_body2prin)

    % magnetic field vector resolved in principal frame [T]
    B = A*B_eci;
    
    % magnetic field vector resolved in body frame [T]
    B_body = R_body2prin'*B;
    
    % sets up gaussian_random_sample MATLAB function for use in Simulink
    coder.extrinsic('gaussian_random_sample');
    
    % initializes noise to avoid Simulink errors
    n = zeros(3,1);
    
    % samples noise [T]
    n = gaussian_random_sample(zeros(3,1),Sigma_mag);
    
    % simulates magnetometer measurement [T]
    Bm_body = B_body+b_mag+n;

end