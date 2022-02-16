%==========================================================================
%
% set_deputy2  Returns the deputy 2 spacecraft parameters.
%
%   deputy2 = set_deputy2
%
% Author: Tamas Kis
% Last Update: 2021-10-18
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   chief       - (struct) chief spacecraft parameters (see 
%                   --> see "initialization/initialize_chief" for full
%                       definition
%
% -------
% OUTPUT:
% -------
%   deputy2     - (struct) deputy 2 spacecraft parameters
%                   --> see "initialization/initialize_deputy" for full
%                       definition
%
%==========================================================================
function deputy2 = set_deputy2(chief)

    % -----------
    % Definition.
    % -----------
    
    % properties
    m = 6;          % mass [kg]
    %AD = 0.09;  	% atmospheric drag reference area [m^2]
    AD = 0.045;  	% atmospheric drag reference area [m^2]
    CD = 2.5;     	% drag coefficient [-]
    Asrp = 0.07; 	% solar radiation pressure reference area [m^2]
    Csrp = 1.32;    % solar radiation pressure coefficient [-]

    % initial relative state (relative orbital elements)
    ada = 0;        % relative semi-major axis [m]
    adlam = 10000;  % relative mean longitude [m]
    adex = 0;       % x-component of relative eccentricity vector [m]
    adey = 1000;    % y-component of relative eccentricity vector [m]
    adix = 0;       % x-component of relative inclination vector [m]
    adiy = 1000;    % y-component of relative inclination vector [m]

    % ---------------
    % Initialization.
    % ---------------
    
    % packages spacecraft properties and relative orbital element state
    properties = [m,AD,CD,Asrp,Csrp];
    ROE = [ada;adlam;adex;adey;adix;adiy];

    % initializes "deputy2" structure
    deputy2 = initialize_deputy(chief,properties,ROE);

end