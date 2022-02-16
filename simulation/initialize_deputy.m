%==========================================================================
%
% initialize_deputy  Initializes the deputy spacecraft's parameters.
%
%   deputy = initialize_deputy(chief,properties,ROE)
%
% Author: Tamas Kis
% Last Update: 2021-08-16
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   chief       - (struct) chief spacecraft parameters (see 
%                 "initialize_chief" for full definition)
%   properties  - (5×1 double) deputy spacecraft properties
%                   1. m    - mass [kg]
%                   2. Ad   - atmospheric drag reference area [m^2]
%                   3. CD   - drag coefficient [-]
%                   4. Asrp - solar radiation pressure reference area [m^2]
%                   5. CR   - coefficient of reflectivity [-]
%   ROE         - (6×1 double) initial relative orbital elements
%                   1. ada      - relative semi-major axis [m]
%                   2. adlam    - relative mean longitude [m]
%                   3. adex     - x-component of relative eccentricity
%                                 vector [m]
%                   4. adey     - y-component of relative eccentricity 
%                                 vector [m]
%                   5. adix     - x-component of relative inclination
%                                 vector [m]
%                   6. adiy     - y-component of relative inclination
%                                 vector [m]
%
% -------
% OUTPUT:
% -------
%   deputy - (struct) deputy spacecraft parameters
%    	• m  	- (1×1 double) mass [kg]
%       • AD   	- (1×1 double) atmospheric drag reference area [m^2]
%    	• CD   	- (1×1 double) drag coefficient [-]
%       • B    	- (1×1 double) ballistic coefficient [-]
%    	• dB  	- (1×1 double) differential ballistic coefficient [-]
%      	• Asrp 	- (1×1 double) solar radiation pressure ref. area [m^2]
%    	• Csrp 	- (1×1 double) solar radiation pressure coefficient [-]
%       • ECI  	- (6×1 double) initial ECI state
%                	1-3. r - position resolved in ECI frame [m]
%                	4-6. v - inertial velocity resolved in ECI frame [m/s]
%       • RTN	- (6×1 double) initial relative RTN state
%                	1-3. dr - rel. position resolved in chief RTN frame [m]
%                	4-6. dv - rel. vel. resolved in chief RTN frame [m/s]
%    	• OE 	- (6×1 double) initial Keplerian orbital element state
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
%               	6. M_m  - mean mean anomaly [rad]
%     	• ROE  	- (6×1 double) initial relative orbital element state
%                   1. ada   - relative semi-major axis [m]
%                   2. adlam - relative mean longitude [m]
%                   3. adex  - x-comp of relative eccentricity vector [m]
%                   4. adey  - y-comp of relative eccentricity vector [m]
%                   3. adix  - x-comp of relative inclination vector [m]
%                   4. adiy  - y-comp of relative inclination vector [m]
%
%==========================================================================
function deputy = initialize_deputy(chief,properties,ROE)

    % extracts spacecraft properties
    m = properties(1);      % mass [kg]
    AD = properties(2);     % atmospheric drag reference area [m^2]
    CD = properties(3);     % drag coefficient [-]
    Asrp = properties(4);   % solar radiation pressure reference area [m^2]
    Csrp = properties(5);   % solar radiation pressure coefficient [-]
    
    % deputy's ballistic coefficient [m^2/kg]
    B = CD*AD/m;
    
    % chief's ballistic coefficient [m^2/kg]
    Bc = chief.B;
    
    % differential ballistic coefficient
    dB = (B-Bc)/Bc;
    
    % initial osculating orbital elements [m][rad]
    %OE = roe2kep(chief.OE,ROE);
    OE = roe2oe(chief.OE,ROE(1),ROE(2),ROE(3),ROE(4),ROE(5),ROE(6));
    
    % initial mean orbital elements [m][rad]
    %MOE = osc2mean(OE);
    MOE = OE;
    
    % Earth gravitational parameter [m^3/s^2]
    %mu = 398600.4415e9;
    
    % initial ECI state [m][m/s]
    %ECI = kep2cart(OE,mu);
    M = OE(6);
    e = OE(2);
    [r_eci,v_eci] = oe2eci(OE(1),OE(2),OE(3),OE(4),OE(5),E2nu(M2E(M,e),e));
    
    ECI = [r_eci;v_eci];

    % initial relative RTN state
    RTN = eci2rtn(chief.ECI,ECI);
    
    % stores parameters in structure
    deputy.M = M;
    deputy.AD = AD;
    deputy.CD = CD;
    deputy.B = B;
    deputy.dB = dB;
    deputy.Asrp = Asrp;
    deputy.CR = Csrp;
    deputy.ECI = ECI;
    deputy.RTN = RTN;
    deputy.OE = OE;
    deputy.MOE = MOE;
    deputy.ROE = ROE;
    
end