%==========================================================================
%
% initialize_chief  Initializes the chief spacecraft's parameters.
%
%   chief = initialize_chief(properties,oe)
%
% Author: Tamas Kis
% Last Update: 2021-10-16
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   properties  - (5×1 double) spacecraft properties
%                   1. m    - mass [kg]
%                   2. Ad   - atmospheric drag reference area [m^2]
%                   3. CD   - drag coefficient [-]
%                   4. Asrp - solar radiation pressure reference area [m^2]
%                   5. CR   - coefficient of reflectivity [-]
%   OE          - (6×1 double) initial osculating orbital elements
%                   1. a  - semi-major axis [m]
%                   2. e  - eccentricity [-]
%                   3. i  - inclination [rad]
%                   4. Om - right ascension of the ascending node [rad]
%                   5. w  - argument of periapsis [rad]
%                   6. M  - mean anomaly [rad]
%
% -------
% OUTPUT:
% -------
%   chief   - (struct) chief spacecraft parameters
%       • m  	- (1×1 double) mass [kg]
%    	• AD   	- (1×1 double) atmospheric drag reference area [m^2]
%    	• CD  	- (1×1 double) drag coefficient [-]
%    	• B   	- (1×1 double) ballistic coefficient [m^2/kg]
%    	• Asrp  - (1×1 double) solar radiation pressure reference area 
%                 [m^2]
%     	• Csrp 	- (1×1 double) solar radiation pressure coefficient [-]
%    	• ECI   - (6×1 double) initial ECI state
%                   1-3. r0_eci - position resolved in ECI frame [m]
%                   4-6. v0_eci - inertial velocity resolved in ECI frame
%                                 [m/s]
%    	• OE  	- (6×1 double) initial Keplerian orbital element state
%                   1. a  - semi-major axis [m]
%                   2. e  - eccentricity [-]
%                 	3. i  - inclination [rad]
%                   4. Om - right ascension of the ascending node [rad]
%                   5. w  - argument of periapsis [rad]
%                   6. M  - mean anomaly [rad]
%    	• MOE 	- (6×1 double) initial mean Keplerian orbital element state
%                	1. a_m  - mean semi-major axis [m]
%                   2. e_m  - mean eccentricity [-]
%                 	3. i_m  - mean inclination [rad]
%                   4. Om_m - mean right ascension of the asc. node [rad]
%                   5. w_m 	- mean argument of periapsis [rad]
%               	6. M    - mean mean anomaly [rad]
%
%==========================================================================
function chief = initialize_chief(properties,OE)

    % extracts spacecraft properties
    m = properties(1);      % mass [kg]
    AD = properties(2);     % atmospheric drag reference area [m^2]
    CD = properties(3);     % drag coefficient [-]
    Asrp = properties(4);   % solar radiation pressure reference area [m^2]
    CR = properties(5);   % solar radiation pressure coefficient [-]
    
    % ballistic coefficient [m^2/kg]
    B = CD*AD/m;
    
    % Earth gravitational parameter [m^3/s^2]
    mu = 398600.4415e9;
    
    % initial ECI state [m][m/s]
    %ECI = kep2cart(OE,mu);
    M = OE(6);
    e = OE(2);
    [r_eci,v_eci] = oe2eci(OE(1),OE(2),OE(3),OE(4),OE(5),E2nu(M2E(M,e),e));
    ECI = [r_eci;v_eci];
    
    % initial mean orbital element state [m][rad]
    %MOE = osc2mean(OE);
    MOE = OE;

    % stores parameters in structure
    chief.M = m;
    chief.AD = AD;
    chief.CD = CD;
    chief.B = B;
    chief.Asrp = Asrp;
    chief.CR = CR;
    chief.ECI = ECI;
    chief.OE = OE;
    chief.MOE = MOE;

end