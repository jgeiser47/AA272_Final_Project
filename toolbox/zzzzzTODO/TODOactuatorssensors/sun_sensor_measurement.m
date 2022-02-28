%==========================================================================
%
% sun_sensor_measurement  Simulates Sun sensor measurements.
%
%   sun_sensor_measurement(b_sun,Sigma_sun,S_hat,R_body2prin)
%
% NOTE: Assumes a multitude of Sun sensors provide Sun-pointing information 
% from incident rays. This measurement function simplifies that behavior by
% assuming these separate sensor outputs are compiled into a measurement of 
% the Sun's unit direction w.r.t. satellite body/payload/sensor frame. This
% way, we don't have to model multiple Sun sensors and how the on-board
% computer combines their respective 2-axis signals. 
%
% Copyright (c) 2021 Luke Neise, Tamas Kis
% Last Update: 2021-06-01
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   S_hat       (3x1) Sun direction (unit vector) w.r.t. satellite resolved
%                     in principal frame
%   b_sun       (3x1) bias
%   Sigma_sun   (3x3) noise covariance
%   R_body2prin (3x3) rotation matrix (body --> principal)
%
% --------
% OUTPUTS:
% --------
%   Sm_body     (3x1) Sun direction measurement resolved in body frame
%
%==========================================================================
function Sm_body = sun_sensor_measurement(S_hat,b_sun,Sigma_sun,...
    R_body2prin)
    
    % resolves S_hat in body frame
    S_hat_body = R_body2prin'*S_hat;

    % sets up gaussian_random_sample MATLAB function for use in Simulink
    coder.extrinsic('gaussian_random_sample');
    
    % initializes noise to avoid Simulink errors
    n = zeros(3,1);
    
    % samples noise
    n = gaussian_random_sample(zeros(3,1),Sigma_sun);

    % simulates Sun sensor measurement
    Sm_body = S_hat_body+b_sun+n;

end