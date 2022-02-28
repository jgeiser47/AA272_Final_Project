% ecef2enu  Determines the position of a satellite relative to a ground
% station in the East, North, Up (ENU) coordinate frame of the ground
% station, given the geodetic coordinates of the ground station and the
% position vector of the satellite in the Earth-Centered, Earth-Fixed
% (ECEF) coordinate frame.
%
%   s_ENU = ecef2enu(r_XYZ,latG,lonG,h) returns the position "s_ENU" of the
%   satellite relative to the ground station in the East, North, Up (ENU)
%   coordinate frame of the ground station, given the geodetic latitude
%   "latG", geodetic longitude "lonG", and geodetic altitude "hG" of the
%   ground station, as well as the position "r_XYZ" of the satellite in the
%   Earth-Centered, Earth-Fixed (ECEF) coordinate frame.
%
% Copyright (c) 2021 Tamas Kis



%% FUNCTION

% INPUT: r_XYZ - position in Earth-Centered, Earth-Fixed (ECEF) coordinate
%                frame
%        latG - geodetic latitude of the ground station [deg]
%        lonG - geodetic longitude of the ground station [deg]
%        hG - geodetic altitude of the ground station [km]
% OUTPUT: s_ENU - position of satellite relative to ground station in the
%                 ground station's East, North, Up (ENU) coordinate frame
function s_ENU = ecef2enu(r_XYZ,latG,lonG,hG)

    % ECEF position of ground station
    rG_XYZ = geod2ecef(latG,lonG,hG);
    
    % precomputes trigonometric functions for speed
    so = sind(latG);
    co = cosd(latG);
    sl = sind(lonG);
    cl = cosd(lonG);
    
    % defines rotation matrix
    R_XYZ_to_ENU = [-sl      cl       0;
                    -so*cl   -so*sl   co;
                    co*cl    co*sl    so];
    
    % position of satellite w.r.t. ground station in ENU coordinate frame
    s_ENU = R_XYZ_to_ENU*(r_XYZ-rG_XYZ);
    
end