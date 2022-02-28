% enu2sae  Determines the range, elevation, and azimuth of a satellite 
% relative to a ground station given the position of the satellite in the
% East, North, Up (ENU) coordinate frame of the ground station.
%
%   [s,A,E] = enu2sae(s_ENU) returns the range "s", azimuth "A", and
%   elevation "E" of a satellite relative to a ground station given the
%   position "s_ENU" of the satellite in the East, North, Up (ENU) 
%   coordinate frame of the ground station.
%
% Copyright (c) 2021 Tamas Kis



%% FUNCTION

% INPUT: s_ENU - position in ground station's East, North, Up (ENU)
%                coordinate frame
% OUTPUT: s - range [km]
%         A - azimuth [deg]
%         E - elevation [deg]
function [s,A,E] = enu2sae(s_ENU)

    % range
    s = norm(s_ENU);
    
    % azimuth
    A = mod(atan2d(s_ENU(1),s_ENU(2)),360);
    
    % elevation
    E = atand(s_ENU(3)/sqrt(s_ENU(1)^2+s_ENU(2)^2));
    
end