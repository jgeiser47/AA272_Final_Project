%==========================================================================
%
% nonlinear_measurement  Nonlinear measurement equation.
% 
% Returns the expected measurement (i.e. the modeled observation vector) 
% using the nonlinear measurement equation.
%
%   z = nonlinear_measurement(x,R_body2prin,B_eci,S_hat_eci,w)
%
% Copyright (c) 2021 Luke Neise, Tamas Kis
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
%--------
% INPUTS:
%--------
%   x           (7x1) prior state vector
%   R_body2prin (3x3) rotation matrix (body frame --> principal frame)
%   B_eci       (3x1) [T] ground truth magnetic field vector resolved in 
%                         ECI frame
%   S_hat_eci   (3x1) Sun direction (unit vector) w.r.t. satellite resolved
%                     in ECI frame
%
%---------
% OUTPUTS:
%---------
%   z           (9x1) expected measurement (modeled observation vector)
%
%==========================================================================
function z = nonlinear_measurement(x,R_body2prin,B_eci,S_hat_eci)

    % unpacks state vector
    q = x(1:4);
    w = x(5:7);
    
    % attitude matrix (ECI --> principal)
    A = quat2dcm(q);

    % expected measurements
    z_mag = R_body2prin'*A*B_eci;
    z_sun = R_body2prin'*A*S_hat_eci;
    z_gyro = R_body2prin'*w;

    % packages modeled observation vector
    z = [z_mag;z_sun;z_gyro];

end