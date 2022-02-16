%==========================================================================
%
% geoc2ecef  ECEF position vector from the geocentric coordinates.
%
%   r_ecef = geoc2ecef(lat_gc,lon_gc,h_gc)
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
%   lat_gc	- (1×1 double) geocentric latitude [deg]
%   lon_gc	- (1×1 double) geocentric longitude [deg]
%   h_gc 	- (1×1 double) geocentric altitude [m]
%   R_sub   - (OPTIONAL) (1×1 double) local Earth radius at the 
%             subsatellite point (default to the mean equatorial radius if 
%             not specified) [m]
%
% -------
% OUTPUT:
% -------
%   r_ecef  - (3×1 double) position resolved in ECEF frame [m]
%
%==========================================================================
function r_ecef = geoc2ecef(lat_gc,lon_gc,h_gc)
    
    % defaults R_sub to the mean equatorial radius [m]
    if nargin == 3
        R_sub = 6378136.3;
    end
    
    % position resolved in ECEF frame [m]
    r_ecef = [(R_sub+h_gc)*cosd(lat_gc)*cosd(lon_gc);
              (R_sub+h_gc)*cosd(lat_gc)*sind(lon_gc);
              (R_sub+h_gc)*sind(lat_gc)];
      
end