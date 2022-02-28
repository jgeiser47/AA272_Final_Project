%==========================================================================
%
% sun_direction  Direction of the Sun with respect to the satellite.
%
%   S_hat = sun_direction(r,r_sun,A)
%
% Author: Tamas Kis
% Last Update: 2021-09-23
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   r       - (3×1 double) satellite position resolved in ECI frame [m]
%   r_sun   - (3×1 double) Sun position resolved in ECI frame [m]
%   A       - (3×3 double) attitude matrix (ECI --> principal) [-]
%
% -------
% OUTPUT:
% -------
%   S_hat   - (3×1 double) direction of Sun (unit vector pointing from 
%             satellite towards Sun) resolved in principal frame
%
%==========================================================================
function S_hat = sun_direction(r,r_sun,A)
    
    % position of Sun w.r.t. satellite, resolved in ECI frame [m]
    S_eci = r_sun-r;
    
    % position of Sun w.r.t. satellite, resolved in principal frame [m]
    S = A*S_eci;
    
    % unit vector corresponding to "S"
    S_hat = S/norm(S);
    
end