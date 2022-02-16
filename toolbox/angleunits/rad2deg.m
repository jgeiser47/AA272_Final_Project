%==========================================================================
%
% rad2deg  Radians to degrees.
%
%   theta_deg = rad2deg(theta_deg)
%
% See also deg2rad, deg2sec, rad2sec, sec2deg, sec2rad.
%
% Author: Tamas Kis
% Last Update: 2022-02-01
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   theta_rad   - (1×1 double) angle in radians [rad]
%
% -------
% OUTPUT:
% -------
%   theta_deg   - (1×1 double) angle in degrees [°]
%
%==========================================================================
function theta_deg = rad2deg(theta_rad)
    theta_deg = (180/pi)*theta_rad;
end