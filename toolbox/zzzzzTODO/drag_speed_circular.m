%==========================================================================
%
% drag_speed_circular  Speed of a spacecraft relative to the atmosphere in
% a near-circular orbit.
%
%   v_atm = drag_speed_circular(a,i,u)
%
% Author: Tamas Kis
% Last Update: 2021-11-14
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   a       - (1×1 double) semi-major axis [m]
%   i       - (1×1 double) inclination [rad]
%   u       - (1×1 double) argument of latitude [rad]
%
% -------
% OUTPUT:
% -------
%   v_atm   - (6×1 double) speed of spacecraft relative to atmosphere [m/s]
%
%==========================================================================
function v_atm = drag_speed_circular(a,i,u)

    % Earth gravitational parameter [m^3/s^2] and rotational speed [rad/s]
    mu_earth = 398600.4415e9;
    w_earth = 0.0000729211585530;
    
    % mean motion [rad/s]
    n = sqrt(mu_earth/a^3);
    
    % terms to simplify calculation of speed relative to atmosphere
    A = 0.75;
    B = (n/w_earth)^2;
    C = -(2*n*cos(i)/w_earth);
    D = (cos(2*i)+cos(2*u))/4;
    E = -((cos(2*i-2*u)+cos(2*i+2*u))/8);
    
    % speed of spacecraft relative to the atmosphere [m/s]
    v_atm = a*w_earth*sqrt(A+B+C+D+E);

end