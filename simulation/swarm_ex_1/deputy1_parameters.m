%==========================================================================
%
% deputy1_parameters  Returns the deputy 1 spacecraft parameters.
%
%   deputy1 = deputy1_parameters
%
% Author: Tamas Kis
% Last Update: 2022-02-25
%
%--------------------------------------------------------------------------
%
% -------
% OUTPUT:
% -------
%   deputy1 - (1×1 struct) deputy 1 satellite parameters
%       • m  	- (1×1 double) mass [kg]
%    	• AD   	- (1×1 double) atmospheric drag reference area [m^2]
%    	• CD  	- (1×1 double) drag coefficient [-]
%    	• B   	- (1×1 double) ballistic coefficient [m^2/kg]
%    	• dB  	- (1×1 double) differential ballistic coefficient [-]
%      	• Asrp 	- (1×1 double) solar radiation pressure ref. area [m^2]
%    	• CR 	- (1×1 double) coefficient of reflectivity [-]
%       • ECI  	- (6×1 double) initial ECI state
%                	1-3. r_eci - position resolved in ECI frame [m]
%                	4-6. v_eci - inertial vel. resolved in ECI frame [m/s]
%       • RTN	- (6×1 double) initial relative RTN state
%                	1-3. dr - rel. position resolved in chief RTN frame [m]
%                	4-6. dv - rel. vel. resolved in chief RTN frame [m/s]
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
%     	• ROE     - (1×1 struct) initial relative orbital elements
%           ‣ ada   - (1×1 double) relative semi-major axis [m]
%           ‣ adlam - (1×1 double) relative mean longitude [m]
%           ‣ adex  - (1×1 double) x-comp. of rel. eccentricity vector [m]
%           ‣ adey  - (1×1 double) y-comp. of rel. eccentricity vector [m]
%           ‣ adix  - (1×1 double) x-comp. of rel. inclination vector [m]
%           ‣ adiy  - (1×1 double) y-comp. of rel. inclination vector [m]
%
%==========================================================================
function deputy1 = deputy1_parameters
    
    % load chief satellite parameters locally
    chief = chief_parameters;

    % properties
    props.m = 6;        % mass [kg]
    props.AD = 0.09;    % atmospheric drag reference area [m^2]
    props.CD = 2.5;     % drag coefficient [-]
    props.Asrp = 0.07;  % solar radiation pressure reference area [m^2]
    props.CR = 1.32;    % coefficient of reflectivity [-]

    % initial relative state (relative orbital elements)
    ROE.ada = 0;        % relative semi-major axis [m]
    ROE.adlam = 1000;   % relative mean longitude [m]
    ROE.adex = 0;       % x-component of relative eccentricity vector [m]
    ROE.adey = 1000;    % y-component of relative eccentricity vector [m]
    ROE.adix = 0;       % x-component of relative inclination vector [m]
    ROE.adiy = 1000;    % y-component of relative inclination vector [m]

    % initializes "deputy1" structure
    deputy1 = initialize_deputy(chief,props,ROE);

end