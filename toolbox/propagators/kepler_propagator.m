%==========================================================================
%
% kepler_propagator  Simulates a Keplerian orbit.
%
%   [t,nu,M,E,r_eci,v_eci,r_ecef,v_ecef,lat,lon,h,GMST] =...
%       kepler_propagator(a,e,i,Om,w,nu0,UTC0_cal,t_length,dt,mu)
%
% Author: Tamas Kis
% Last Update: 2021-12-12
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   a           - (1×1 double) semi-major axis [m]
%   e           - (1×1 double) eccentricity [-]
%   i           - (1×1 double) inclination [rad]
%   Om          - (1×1 double) right ascension of the ascending node [rad]
%   w           - (1×1 double) argument of periapsis [rad]
%   nu0         - (1×1 double) true anomaly at epoch [rad]
%   UTC0_cal    - (1×6 double) epoch (UTC time) [YYYY,MM,DD,mm,hh,ss]
%   t_length    - (OPTIONAL) (1×1 double) simulation length (defaults to
%                 one orbital period) [s]
%   dt          - (OPTIONAL) (1×1 double) time step (defaults to 
%                 t_length/1000) [s]
%   mu          - (OPTIONAL) (1×1 double) gravitational parameter [m^3/s^2]
%
% -------
% OUTPUT:
% -------
%   t           - (N×1 double) simulation time vector [s]
%   nu          - (N×1 double) true anomaly [rad]
%   M           - (N×1 double) mean anomaly [rad]
%   E           - (N×1 double) eccentric anomaly [rad]
%   r_eci       - (N×3 double) position resolved in ECI frame [m]
%   v_eci       - (N×3 double) inertial vel. resolved in ECI frame [m/s]
%   r_ecef      - (N×3 double) position resolved in ECEF frame [m]
%   v_ecef      - (N×3 double) ECEF velocity resolved in ECEF frame [m/s]
%   lat         - (N×1 double) geodetic latitude [deg]
%   lon         - (N×1 double) geodetic longitude [deg]
%   h           - (N×1 double) geodetic altitude [m]
%   GMST        - (N×1 double) Greenwich mean sidereal time [rad]
%
% -----
% NOTE:
% -----
%   --> time vector (t) is measured with respect to the epoch (i.e. t = 1 s
%       means one second after epoch)
%   --> N = length of time vector
%   
%==========================================================================
function [t,nu,M,E,r_eci,v_eci,r_ecef,v_ecef,lat,lon,h,GMST] =...
    kepler_propagator(a,e,i,Om,w,nu0,UTC0_cal,t_length,dt,mu)
    
    % ----------------------------------------------------
    % Sets unspecified parameters to their default values.
    % ----------------------------------------------------

    % sets gravitational parameter (defaults to Earth's) [m^3/s^2]
    if nargin < 10
        mu = 398600.4415e9;
    end

    % sets simulation length (defaults to one orbital period) [s]
    if (nargin < 8) || isempty(t_length)
        t_length = a2T(a,mu);
    end

    % sets time step (defaults to t_length/1000) [s]
    if (nargin < 9) || isempty(dt)
        dt = round(t_length/1000);
    end

    % ------------------
    % Orbit propagation.
    % ------------------
    
    % simulation time vector
    t = (0:dt:t_length)';

    % converts epoch (UTC time) to MJD
    UTC0_MJD = cal2mjd(UTC0_cal);
    
    % difference between UTC and UT1 time (DUT1 = UT1 - UTC) [s]
    DUT1 = getDUT1(UTC0_cal);
    
    % gets epoch in UT1 time
    UT10_MJD = UTCtoUT1(UTC0_MJD,DUT1);
    
    % mean motion
    n = a2n(a,mu);
    
    % mean anomaly at epoch
    M0 = nu2M(nu0,e);
    
    % preallocates arrays
    nu = zeros(size(t));
    M = zeros(size(t));
    E = zeros(size(t));
    r_eci = zeros(3,length(t));
    v_eci = zeros(3,length(t));
    r_ecef = zeros(3,length(t));
    v_ecef = zeros(3,length(t));
    lat = zeros(size(t));
    lon = zeros(size(t));
    h = zeros(size(t));
    GMST = zeros(size(t));
    
    % orbit propagation
    for j = 1:length(t)
        
        % mean anomaly from time
        M(j) = t2M(t(j),t(1),M0,n);
        
        % eccentric anomaly from mean anomaly
        E(j) = M2E(M(j),e);
        
        % true anomaly from eccentric anomaly
        nu(j) = E2nu(E(j),e);
        
        % ECI position [m] and velocity [m/s] from orbital elements
        [r_eci(:,j),v_eci(:,j)] = oe2eci(a,e,i,Om,w,nu(j));
        
        % obtains current UT1 time in MJD
        UT1_MJD = t(j)/(24*3600)+UT10_MJD;
        
        % Greenwich mean sidereal time [rad] from UT1 time
        GMST(j) = UT1toGMST(UT1_MJD);
        
        % ECEF position [m] and velocity [m/s]
        [r_ecef(:,j),v_ecef(:,j)] = eci2ecef(r_eci(:,j),v_eci(:,j),...
            GMST(j));
        
        % geodetic coordinates from ECEF position
        [lat(j),lon(j),h(j)] = ecef2geod(r_ecef(:,j));
        
    end

end