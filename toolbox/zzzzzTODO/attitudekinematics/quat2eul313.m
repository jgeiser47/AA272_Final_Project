%==========================================================================
%
% quat2eul313  Quaternion to 3-1-3 Euler angles.
%
%   [phi,theta,psi] = quat2eul313(q)
%
% Author: Tamas Kis
% Last Update: 2021-08-26
%
% REFERENCES:
%   [1] https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Quaternion_%E2%86%92_Euler_angles_(z-x-z_extrinsic)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   q       - (4×1 double) quaternion [-]
%
% -------
% OUTPUT:
% -------
%   phi     - (1×1 double) 3-1-3 Euler angle [rad]
%   theta   - (1×1 double) 3-1-3 Euler angle [rad]
%   psi     - (1×1 double) 3-1-3 Euler angle [rad]
%
%==========================================================================
function [phi,theta,psi] = quat2eul313(q)

    % unpacks quaternion
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    
    % 3-1-3 Euler angles [rad]
    phi = mod(atan2(q1*q3+q2*q4,q1*q4-q2*q3),2*pi);
    theta = acos(-q1^2-q2^2+q3^2+q4^2);
    psi = mod(atan2(q1*q3-q2*q4,q2*q3+q1*q4),2*pi);
    
end