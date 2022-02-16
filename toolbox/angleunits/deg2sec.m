%==========================================================================
%
% deg2sec  Degrees to arcseconds.
%
%   theta_arcsec = sec2deg(theta_deg)
%
% See also deg2rad, rad2deg, rad2sec, sec2deg, sec2rad.
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
%   theta_deg       - (1×1 double) angle in degrees [°]
%
% -------
% OUTPUT:
% -------
%   theta_arcsec    - (1×1 double) angle in arcseconds ['']
%
%==========================================================================
function theta_arcsec = deg2sec(theta_deg)
    theta_arcsec = 3600*theta_deg;
end