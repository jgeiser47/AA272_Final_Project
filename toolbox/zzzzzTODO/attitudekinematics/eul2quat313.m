%==========================================================================
%
% eul2quat313  3-1-3 Euler angles to quaternion.
%
%   q = eul2quat313(phi,theta,psi)
%
% Author: Tamas Kis
% Last Update: 2021-09-09
%
% REFERENCES:
%   [1] Wertz, "Spacecraft Attitude Determination and Control", Eq. (E-10)
%       (p. 762)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   phi   	- (1×1 double) 3-1-3 Euler angle [rad]
%   theta 	- (1×1 double) 3-1-3 Euler angle [rad]
%   psi     - (1×1 double) 3-1-3 Euler angle [rad]
%
% -------
% OUTPUT:
% -------
%   q       - (4×1 double) quaternion [-]
%
%==========================================================================
function q = eul2quat313(phi,theta,psi)

    % calculates quaternion elements
    q1 = sin(theta/2)*cos((phi-psi)/2);
    q2 = sin(theta/2)*sin((phi-psi)/2);
    q3 = cos(theta/2)*sin((phi+psi)/2);
    q4 = cos(theta/2)*cos((phi+psi)/2);
    
    % packages quaternion
    q = [q1;q2;q3;q4];
 
end