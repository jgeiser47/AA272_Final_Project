%==========================================================================
%
% r2ad  ECI position to geocentric right ascension and declination.
%
%   [alpha,delta] = r2ad(r_eci)
%
% Author: Tamas Kis
% Last Update: 2022-01-06
%
% REFERENCES:
%   [1] Vallado, "Fundamentals of Astrodynamics and Applications", 4th Ed.
%       (pp. 256-259)
%   [2] Montenbruck and Gill, "Satellite Orbits" (p. 25)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   r_eci   - (3×1 double) position in resolved ECI frame [any unit]
%
% -------
% OUTPUT:
% -------
%   alpha   - (1×1 double) right ascension [rad]
%   delta   - (1×1 double) declination [rad]
%
%==========================================================================
function [alpha,delta] = r2ad(r_eci)
    
    % extracts components of position vector
    rI = r_eci(1);
    rJ = r_eci(2);
    rK = r_eci(3);
    
    % right ascension [rad]
    alpha = iatan2(rJ,rI);
    
    % declination [rad]
    delta = iatan2(rK,sqrt(rI^2+rJ^2));
    
end