%==========================================================================
%
% ecef2geod  Geodetic coordinates from the ECEF position vector.
%
%   [lat,lon,h] = ecef2geod(r_ecef)
%   [lat,lon,h] = ecef2geod(r_ecef,TOL)
%   [lat,lon,h] = ecef2geod(r_ecef,[],imax)
%   [lat,lon,h] = ecef2geod(r_ecef,TOL,imax)
%
% Author: Tamas Kis
% Last Update: 2021-10-19
%
% REFERENCES:
%   [1] Vallado, "Fundamentals of Astrodynamics and Applications", 4th Ed.
%       (pp. 138, 169-172)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   r_ecef  - (3×1 double) position resolved in ECEF frame [m]
%   TOL     - (1×1 double) (OPTIONAL) tolerance (defaults to 1e-12)
%   imax    - (1×1 double) (OPTIONAL) maximum number of iterations 
%             (defaults to 1e6)
%
% -------
% OUTPUT:
% -------
%   lat     - (1×1 double) geodetic latitude [deg]
%   lon    	- (1×1 double) geodetic longitude [deg]
%   h      	- (1×1 double) geodetic altitude [m]
%
%==========================================================================
function [lat,lon,h] = ecef2geod(r_ecef,TOL,imax)
    
    % defaults tolerance and maximum number of iterations if not input
    if (nargin < 2) || isempty(TOL)
        TOL = 1e-12;
    end
    if (nargin < 3) || isempty(imax)
        imax = 1e6;
    end
        
    % mean equatorial radius [m] and eccentricity [-] of the Earth
    R_earth = 6378136.3;
    e_earth = 0.081819221456;
    
    % extracts X, Y, and Z coordinates of ECEF position vector [m]
    rX = r_ecef(1);
    rY = r_ecef(2);
    rZ = r_ecef(3);
    
    % magnitude of position vector [m]
    r_mag = norm(r_ecef);
    
    % scalar projection of position vector onto equatorial plane [m]
    rXY = sqrt(rX^2+rY^2);
    
    % geodetic longitude [deg]
    %lon = atan2d(rY,rX);
    lon = iatan2d(rY,rX);
    
    % initial guess for geodetic latitude [deg]
    lat_old = asind(rZ/r_mag);
    
    % initializes the error so the loop will be entered
    err = 2*TOL;
    
    % initalizes radius of curvature in the meridian (needed for Simulink)
    N_phi = 0;
    
    % fixed-point iteration to solve for geodetic latitude
    i = 1;
    while (err > TOL) && (i < imax)

        % radius of curvature in the meridian [m]
        N_phi = R_earth/sqrt(1-e_earth^2*sind(lat_old)^2);
        
        % updates estimate of geodetic latitude [deg]
        %lat_new = atan2d(rZ+N_phi*e_earth^2*sind(lat_old),rXY);
        lat_new = iatan2d(rZ+N_phi*e_earth^2*sind(lat_old),rXY);

        % calculates error
        err = abs(lat_new-lat_old);

        % stores updated estimate of geodetic latitude for next iteration
        lat_old = lat_new;

        % increments loop index
        i = i+1;

    end

    % geodetic latitude [deg]
    lat = lat_new;
    
    % geodetic altitude [m]
    h = rXY/cosd(lat)-N_phi;
      
end