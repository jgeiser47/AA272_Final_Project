%==========================================================================
%
% w_earth_iau06  Angular velocity of the Earth resolved in the ITRF
% (IAU2006/2000, CIO based).
%
%   w_itrf = w_earth_iau06(LOD,W)
%
% See also w_earth_approx.
%
% Author: Tamas Kis
% Last Update: 2022-02-15
%
% REFERENCES:
%   [1] https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html
%   [2] Vallado, "Fundamentals of Astrodynamics and Applications", 4th Ed.
%       (p. 220)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   LOD     - (1×1 double) length of day [ms]
%   W       - (3×3 double) polar-motion matrix (ITRF --> TIRS)
%
% -------
% OUTPUT:
% -------
%   w_itrf  - (3×1 double) Earth angular velocity resolved in ITRF [rad/s]
%
%==========================================================================
function w_itrf = w_earth_iau06(LOD,W)
    
    % convert length of day to seconds
    LOD = 0.001*LOD;

    % rotation rate of the Earth [rad/s]
    w_earth = (7.292115146706979e-5)*(1-(LOD/86400));

    % angular velocity of the Earth resolved in the TIRS frame [rad/s]
    w_TIRS = [0;
              0;
              w_earth];

    % angular velocity of the Earth resolved in the ITRF [rad/s]
    w_itrf = (W.')*w_TIRS;
    
end