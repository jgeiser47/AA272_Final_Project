%==========================================================================
%
% oe2eci  Keplerian orbital elements to ECI position and velocity.
%
%   [r_eci,v_eci] = oe2eci(a,e,i,Om,w,nu)
%
% Author: Tamas Kis
% Last Update: 2021-08-10
%
% REFERENCES:
%   [1] Vallado, "Fundamentals of Astrodynamics and Applications", 4th Ed.,
%       Algorithm 10 (pp. 118-119)
%   [2] D'Amico, "Keplerian Orbital Elements, Coordinate Systems", AA 279A 
%       Lecture 5 Slides, Slide 8
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   a       - (1×1 double) semi-major axis [m]
%   e       - (1×1 double) eccentricity [-]
%   i       - (1×1 double) inclination [rad]
%   Om      - (1×1 double) right ascension of the ascending node [rad]
%   w       - (1×1 double) argument of periapsis [rad]
%   nu      - (1×1 double) true anomaly [rad]
%
% -------
% OUTPUT:
% -------
%   r_eci   - (3×1 double) position resolved in ECI frame [m]
%   v_eci 	- (3×1 double) inertial velocity resolved in ECI frame [m/s]
%
% ---------------------
% HANDLING SINGULARIES:
% ---------------------
%  (1) non-circular + equatorial --> input Om = 0 and w = Pi, where Pi is
%                                    the true longitude of periapsis
%  (2) circular + inclined --> input w = 0 and nu = u, where u is the true
%                              argument of latitude
%  (3) circular + equatorial --> input w = 0, Om = 0, and nu = l, where l 
%                                is the true longitude
%
%==========================================================================
function [r_eci,v_eci] = oe2eci(a,e,i,Om,w,nu)
    
    % Earth gravitational parameter [m^3/s^2]
    mu_earth = 398600.4415e9;
    
    % eccentric anomaly [rad]
    E = mod(2*atan(sqrt((1-e)/(1+e))*tan(nu/2)),2*pi);
    
    % mean motion [rad/s]
    n = sqrt(mu_earth/a^3);
    
    % position [m] and inertial velocity [m/s] resolved in PQW frame
    r_pqw = [a*(cos(E)-e);
             a*sqrt(1-e^2)*sin(E);
             0];
    v_pqw = ((a*n)/(1-e*cos(E)))*[-sin(E);
                                   sqrt(1-e^2)*cos(E);
                                   0];
    
    % precomputes trigonometric functions
    sO = sin(Om);
    cO = cos(Om);
    si = sin(i);
    ci = cos(i);
    sw = sin(w);
    cw = cos(w);
    
    % defines rotation matrix
    R_pqw2eci = [cO*cw-sO*ci*sw   -cO*sw-sO*ci*cw    sO*si;
                 sO*cw+cO*ci*sw   -sO*sw+cO*ci*cw   -cO*si;
                 si*sw             si*cw             ci];
    
    % position [m] and inertial velocity [m/s] resolved in ECI frame
    r_eci = R_pqw2eci*r_pqw;
    v_eci = R_pqw2eci*v_pqw;
    
end