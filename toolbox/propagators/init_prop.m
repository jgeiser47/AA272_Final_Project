%==========================================================================
%
% init_prop  Propagator initialization.
%
%   prop = init_prop(models,perturb,UTC_start,duration,time_step,...
%       integrator,sim_start)
%
% Author: Tamas Kis
% Last Update: 2022-02-15
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   models      - (1×1 struct) environmental models for propagator
%       • density     - (char) specifies density model ('Harris-Priester', 
%                       'Jacchia-Roberts', or 'Exponential')
%       • grav_mu     - (1×1 double) Earth gravitational parameter for
%                       gravity model [m^3/s^2]
%       • grav_N      - (1×1 double) degree/order of gravity model
%       • grav_R      - (1×1 double) Earth mean equatorial radius for
%                       gravity model [m]
%   perturb     - (1×1 struct) turns perturbations on (1) or off (0)
%       • drag  	 - (1×1 double or logical) atmospheric drag
%       • emp_accel	 - (1×1 double or logical) empirical acceleration
%       • moon  	 - (1×1 double or logical) 3rd body gravity of Moon
%       • relativity - (1×1 double or logical) relativistic effects
%       • srp    	 - (1×1 double or logical) solar radiation pressure
%       • sun   	 - (1×1 double or logical) 3rd body gravity of Sun
%   UTC_start   - (1×1 double) initial UTC [y,mo,d,h,m,s]
%   duration    - (1×1 double) simulation duration [hr]
%   time_step   - (1×1 double) time step [s]
%   integrator  - (char) integrator ('RK4' or 'ABM8')
%   sim_start   - (1×1 double) initial simulation time [s]
%
% -------
% OUTPUT:
% -------
%   prop        - (1×1 struct) structure to store all propagator settings
%       • data - (1×1 struct) TODO
%     	• models - (1×1 struct) environmental models for propagator (same
%                  as input)
%      	• perturb - (1×1 struct) turns perturbations on (1) or off (0)
%      	            (same as input)
%     	• t0 - (1×1 struct) initial times
%           ‣ t0       - (1×1 double) simulation time [s]
%           ‣ MJD_GPS0 - (1×1 double) GPS time [MJD]
%           ‣ MJD_TAI0 - (1×1 double) TAI (International Atomic Time) [MJD]
%           ‣ MJD_TT0  - (1×1 double) TT (Terrestrial Time) [MJD]
%           ‣ MJD_UT10 - (1×1 double) UT1 (Universal Time 1) [MJD]
%           ‣ MJD_UTC0 - (1×1 double) UTC (Universal Coordinated Time) 
%                        [MJD]
%
%==========================================================================
function prop = init_prop(models,perturb,UTC_start,duration,time_step,...
    integrator,sim_start)

    % --------------
    % Initial times.
    % --------------

    % defaults initial simulation time to 0 if not input
    if (nargin < 6) || isempty(sim_start)
        sim_start = 0;
    end

    % modified Julian date of initial UTC [MJD]
    MJD_UTC0 = cal2mjd(UTC_start);

    % difference between UT1 and UTC (ΔUT1 = UT1 - UTC) [s]
    DUT1 = get_DUT1(MJD_UTC0);

    % difference between TAI and UTC (ΔAT = TAI - UTC) [s]
    DAT = get_DAT(MJD_UTC0);

    % modified Julian dates of initial times in remaining time scales [MJD]
    MJD_UT10 = utc2ut1(MJD_UTC0,DUT1);
    MJD_TAI0 = utc2tai(MJD_UTC0,DAT);
    MJD_GPS0 = tai2gps(MJD_TAI0);
    MJD_TT0 = tai2tt(MJD_TAI0);

    % packages initial times into structure
    t0.t0 = sim_start;
    t0.MJD_GPS0 = MJD_GPS0;
    t0.MJD_TAI0 = MJD_TAI0;
    t0.MJD_TT0 = MJD_TT0;
    t0.MJD_UT10 = MJD_UT10;
    t0.MJD_UTC0 = MJD_UTC0;

    % --------------------
    % ODE solver settings.
    % --------------------

    % integration method
    ode_solver.method = integrator;

    % simulation start time [s]
    ode_solver.t0 = t0.t0;

    % simulation end time [s]
    ode_solver.tf = ode_solver.t0+(duration*3600);

    % time step [s]
    ode_solver.dt = time_step;

    % -----
    % Data.
    % -----

    % gravity coefficients
    GGM05S = load_GGM05S(models.grav_N);

    % function handles for returning normalized grarvity coefficients
    [data.C,data.S] = gravity_coeffs(GGM05S.C_norm,GGM05S.S_norm);

    % Earth orientation parameters
    data.eop = load_eop(MJD_UTC0,duration);

    % additional data for calculating Earth orientation
    data.XYs_iau06 = load_XYs_iau06;
    
    % space weather data
    data.Kp_ap_F107 = load_Kp_ap_F107(MJD_UTC0,duration);

    % -------------------------------------
    % Package overall propagator structure.
    % -------------------------------------

    prop.data = data;
    prop.models = models;
    prop.ode_solver = ode_solver;
    prop.perturb = perturb;
    prop.t0 = t0;

end