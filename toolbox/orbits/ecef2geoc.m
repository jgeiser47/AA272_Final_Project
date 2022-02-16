%==========================================================================
%
% ecef2geoc  Geocentric coordinates from the ECEF position vector.
%
%   [lat_gc,lon_gc,h_gc] = ecef2geoc(r_ecef)
%
% Author: Tamas Kis
% Last Update: 2021-10-05
%
% REFERENCES:
%   [1] D'Amico, "Earth Rotation and Time Systems", AA 279A Lecture 6
%       Slides (p. 10)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   r_ecef  - (3×1 double) position resolved in ECEF frame [m]
%   R_sub   - (OPTIONAL) (1×1 double) local Earth radius at the 
%             subsatellite point (default to the mean equatorial radius if 
%             not specified) [m]
%
% -------
% OUTPUT:
% -------
%   lat_gc  - (1×1 double) geodetic latitude [deg]
%   lon_gc	- (1×1 double) geodetic longitude [deg]
%   h_gc 	- (1×1 double) geodetic altitude [m]
%
%==========================================================================
function [lat_gc,lon_gc,h_gc] = ecef2geoc(r_ecef,R_sub)

    % defaults R_sub to the mean equatorial radius [m]
    if nargin == 1 
        R_sub = 6378136.3;
    end
    
    % extracts X, Y, and Z coordinates of position vector [km]
    rX = r_ecef(1);
    rY = r_ecef(2);
    rZ = r_ecef(3);
    
    % magnitude of position vector [km]
    r = norm(r_ecef);
    
    % geocentric latitude [deg], longitude [deg], and altitude [km]
    lat_gc = asind(rZ/r);
    lon_gc = atan2d(rY,rX);
    h_gc = r-R_sub;
    
end