%==========================================================================
%
% simulate_orbit  Simulates an orbit.
%
%   data = simulate_orbit(spacecraft,prop)
%
% Author: Tamas Kis
% Last Update: 2021-11-15
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   spacecraft  - (struct) spacecraft parameters, see "settings/
%                 set_chief.m", "settings/set_deputy1.m", etc.
%   prop        - (struct) propagator settings, see "settings/
%                 set_propagator.m"
%
% -------
% OUTPUT:
% -------
%   data        - (struct) structure to store all simulation data
%       • r_eci    - (3×N double) position resolved in ECI frame [m]
%       • v_eci    - (3×N double) inertial vel. resolved in ECI frame [m/s]
%       • r_ecef   - (3×N double) position resolved in ECEF frame [m]
%       • v_ecef   - (3×N double) ECEF vel. resolved in ECEF frame [m/s]
%       • lat      - (N×1 double) geodetic latitude [deg]
%       • lon      - (N×1 double) geodetic longitude [deg]
%       • h        - (N×1 double) geodetic altitude [deg]
%       • rho_nrlm - (N×1 double) atmospheric density from NRLMSISE-00
%                    model [kg/m^3]
%       • rho_hp   - (N×1 double) atmospheric density from Harris-Priester 
%                    model [kg/m^3]
%
% -----
% NOTE:
% -----
%   --> N = number of iterations (i.e. length of time vector)
%
%--------------------------------------------------------------------------
function data = simulate_orbit(sat,prop)

    % --------------------------
    % Setup of orbit simulation.
    % --------------------------
    
    % extracts data needed for orbit propagation
    %eop = prop.data.eop;
    %leap_second = prop.data.leap_second;
    %gravity_coeffs = prop.data.gravity_coeff;
    %Kp_ap_F107 = prop.data.Kp_ap_F107;
    
    % extracts needed parameters TODO
    dt = prop.ode_solver.dt;
    x0 = sat.ECI;
    t0 = prop.ode_solver.t0;
    tf = prop.ode_solver.tf;
    
    % dynamics equation
    f = @(t,x) newton_propagator(t,x,prop,sat);

    % simulates orbit
    if strcmpi(prop.ode_solver.method,'ABM8')
        [t,x] = ABM8(f,[t0,tf],x0,dt,true);
    elseif strcmpi(prop.ode_solver.method,'RK4')
        [t,x] = RK4(f,[t0,tf],x0,dt,true);
    end

    x = x';

    data.t = t;
    
    data.r_eci = x(1:3,:);
    data.v_eci = x(4:6,:);


    data.a_eci = zeros(3,length(t));
    data.r_ecef = zeros(3,length(t));
    data.v_ecef = zeros(3,length(t));
    data.lat = zeros(length(t),1);
    data.lon = zeros(length(t),1);
    data.h = zeros(length(t),1);
    data.rho = zeros(length(t),1);

    data.g_eci = zeros(3,length(t));

    for i = 1:length(t)
        [~,extra] = f(t(i),x(:,i));
        data.a_eci(:,i) = extra.a_eci;
        data.r_ecef(:,i) = extra.r_ecef;
        data.v_ecef(:,i) = extra.v_ecef;
        data.lat(i) = extra.lat;
        data.lon(i) = extra.lon;
        data.h(i) = extra.h;
        data.rho(i) = extra.rho;

        data.g_eci(:,i) = extra.g_eci;
    end
    
end