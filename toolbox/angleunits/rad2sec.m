%==========================================================================
%
% rad2sec  Radians to arcseconds.
%
%   theta_sec = rad2sec(theta_rad)
%
% See also deg2rad, deg2sec, rad2deg, sec2deg, sec2rad.
%
% Author: Tamas Kis
% Last Update: 2022-02-01
%
% REFERENCES:
%   [1] Vallado, "Fundamentals of Astrodynamics and Applications", 4th Ed.
%       (p. 175)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   theta_rad       - (1×1 double) angle in radians [rad]
%
% -------
% OUTPUT:
% -------
%   theta_arcsec    - (1×1 double) angle in arcseconds ['']
%
%==========================================================================
function theta_arcsec = rad2sec(theta_rad)
    theta_arcsec = (648000/pi)*theta_rad;
end