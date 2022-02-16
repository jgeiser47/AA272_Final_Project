%==========================================================================
%
% geod2ecef  ECEF position vector from the geodetic coordinates.
%
%   r_ecef = geod2ecef(lat,lon,h)
%
% Author: Tamas Kis
% Last Update: 2021-11-14
%
% REFERENCES:
%   [1] D'Amico, "Plotting of Orbits and Ground Tracks", AA 279A Lecture 
%       7-8 Slides (p. 6)
%   [2] Vallado, "Fundamentals of Astrodynamics and Applications", 4th Ed.,
%       p. 144
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   lat     - (1×1 double) geodetic latitude [deg]
%   lon    	- (1×1 double) geodetic longitude [deg]
%   h      	- (1×1 double) geodetic altitude [m]
%
% -------
% OUTPUT:
% -------
%   r_ecef  - (3×1 double) position resolved in ECEF frame [m]
%
%==========================================================================
function r_ecef = geod2ecef(lat,lon,h)
    
    % mean equatorial radius [m] and eccentricity [-] of the Earth
    R_earth = 6378136.3;
    e_earth = 0.0818;

    % radius of curvature in the meridian [m]
    N_phi = R_earth/sqrt(1-e_earth^2*sind(lat)^2);
    
    % position resolved in ECEF frame [m]
    r_ecef = [(N_phi+h)*cosd(lat)*cosd(lon);
              (N_phi+h)*cosd(lat)*sind(lon);
              (N_phi*(1-e_earth^2)+h)*sind(lat)];
      
end