%==========================================================================
%
% chief_parameters  Define the chief satellite's parameters.
%
%   chief = chief_parameters
%
% Author: Tamas Kis
% Last Update: 2022-02-25
%
%--------------------------------------------------------------------------
%
% -------
% OUTPUT:
% -------
%   chief   - (1×1 struct) chief satellite parameters
%       • m  	- (1×1 double) mass [kg]
%    	• AD   	- (1×1 double) atmospheric drag reference area [m^2]
%    	• CD  	- (1×1 double) drag coefficient [-]
%    	• B   	- (1×1 double) ballistic coefficient [m^2/kg]
%    	• Asrp  - (1×1 double) solar radiation pressure reference area 
%                 [m^2]
%     	• Csrp 	- (1×1 double) coefficient of reflectivity [-]
%    	• ECI   - (6×1 double) initial ECI state
%                   1-3. r_eci - position resolved in ECI frame [m]
%                   4-6. v_eci - inertial velocity resolved in ECI frame
%                                [m/s]
%    	• OE  	- (1×1 struct) initial Keplerian orbital element state
%           ‣ a  - (1×1 double) semi-major axis [m]
%           ‣ e  - (1×1 double) eccentricity [-]
%           ‣ i  - (1×1 double) inclination [rad]
%           ‣ Om - (1×1 double) right ascension of the ascending node [rad]
%           ‣ w  - (1×1 double) argument of periapsis [rad]
%           ‣ nu - (1×1 double) true anomaly [rad]
%    	• MOE 	- (1×1 struct) initial mean Keplerian orbital element state
%           ‣ a_m  - (1×1 double) mean semi-major axis [m]
%           ‣ e_m  - (1×1 double) mean eccentricity [-]
%           ‣ i_m  - (1×1 double) mean inclination [rad]
%           ‣ Om_m - (1×1 double) mean RAAN [rad]
%           ‣ w_m  - (1×1 double) mean argument of periapsis [rad]
%           ‣ M    - (1×1 double) mean mean anomaly [rad]
%
%==========================================================================
function chief = chief_parameters
    
    % properties
    props.m = 6;            % mass [kg]
    props.AD = 0.01;        % atmospheric drag reference area [m^2]
    props.CD = 2.5;         % drag coefficient [-]
    props.Asrp = 0.07;      % solar radiation pressure reference area [m^2]
    props.CR = 1.32;        % coefficient of reflectivity [-]

    % initial absolute state (osculating orbital elements)
    OE.a = 6853137;         % semi-major axis [m]
    OE.e = 0.001;           % eccentricity [-]
    OE.i = deg2rad(51.6);   % inclination [rad]
    OE.Om = 0;              % right ascension of the ascending node [rad]
    OE.w = 0;               % argument of periapsis [rad]
    OE.nu = pi/2;           % true anomaly [rad]

    % initializes "chief" structure
    chief = initialize_chief(props,OE);

end