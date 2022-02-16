%==========================================================================
%
% chief  Define the chief satellite's parameters.
%
%   chief = set_chief
%
% Author: Tamas Kis
% Last Update: 2022-02-07
%
%--------------------------------------------------------------------------
%
% -------
% OUTPUT:
% -------
%   chief   - (struct) chief spacecraft parameters
%               --> see "initialization/initialize_chief" for full
%                   definition
%
%==========================================================================
function chief = set_chief

    % -----------
    % Definition.
    % -----------
    
    % properties
    m = 6;          % mass [kg]
    AD = 0.01;     	% atmospheric drag reference area [m^2]
    CD = 2.5;     	% drag coefficient [-]
    Asrp = 0.07; 	% solar radiation pressure reference area [m^2]
    Csrp = 1.32;    % solar radiation pressure coefficient [-]

    % chief initial absolute state (osculating orbital elements)
    a = 6853137;    % semi-major axis [m]
    e = 0.001;      % eccentricity [-]
    i = 0.9006;     % inclination [rad]
    Om = 0;         % right ascension of the ascending node [rad]
    w = 0;          % argument of periapsis [rad]
    M = pi/2;       % mean anomaly [rad]

    % ---------------
    % Initialization.
    % ---------------
    
    % packages spacecraft properties and orbital element state
    properties = [m;AD;CD;Asrp;Csrp];
    OE = [a;e;i;Om;w;M];

    % initializes "chief" structure
    chief = initialize_chief(properties,OE);

end