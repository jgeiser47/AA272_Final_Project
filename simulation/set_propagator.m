%==========================================================================
%
% set_propagator  Propagator settings.
%
%   prop = set_propagator
%
% Author: Tamas Kis
% Last Update: 2021-10-18
%
%--------------------------------------------------------------------------
%
% -------
% OUTPUT:
% -------
%   prop    - (struct) structure to store all propagator settings
%       • data   - (struct) see "initialization/initialize_data.m"
%     	• models - (struct) stores model settings
%           ‣ control_frame    - (1×1 double) specifies reference frame for
%                                control actions
%           ‣ density          - (1×1 double) specifies density model
%           ‣ gravity_order    - (1×1 double) order of gravity model
%           ‣ gravity_degree   - (1×1 double) degree of gravity model
%           ‣ integrator       - (1×1 double) specifies integrator
%           ‣ reference_system - (1×1 double) specifies reference system
%                                for Earth orientation
%      	• perturbations - (struct) turns perturbations on (1) or off (0)
%           ‣ drag  	 - (1×1 double) atmospheric drag
%           ‣ emp_accel	 - (1×1 double) empirical acceleration
%           ‣ moon  	 - (1×1 double) 3rd body gravity of Moon
%           ‣ polar    	 - (1×1 double) polar motion of Earth's axis
%           ‣ relativity - (1×1 double) relativistic effects
%           ‣ srp    	 - (1×1 double) solar radiation pressure
%           ‣ sun   	 - (1×1 double) 3rd body gravity of Sun
%     	• time - (struct) see "initialization/initialize_time.m"
%
%==========================================================================
function prop = set_propagator

    % ---------
    % Settings.
    % ---------
    
    % epoch (UTC) [YYYY,MM,DD,hh,mm,ss]
    UTC0_cal = [2017,1,1,0,0,0];
    
    % timing
    t_duration = 48;    % simulation duration [hr]
    dt = 10;            % time step [s]
    %t_duration = 5;    % simulation duration [hr]
    %dt = 1;            % time step [s]
    
    % model settings
	models.control_frame = 0;       % control frame (0 = ECI, 1 = RTN)
    models.density = 0;             % density model (0 = NRLMSISE-00,
                                    % 1 = Harris-Priester)
    models.gravity_degree = 120;  	% degree of gravity model [-]
    models.gravity_order = 120;    	% order of gravity model [-]
    models.integrator = 0;          % integrator (0 = RK4, 1 = RK78)
    models.reference_system = 0;    % reference system (0 = IAU 2006/2000A,
                                    % 1 = EME2000)
    
    % perturbations to apply (1 = yes, 0 = no)
    perturbations.drag = 1;         % atmospheric drag
    perturbations.emp_accel = 0;    % empirical acceleration
	perturbations.moon = 1;         % 3rd body gravity of Moon
    perturbations.polar = 1;      	% polar motion of Earth's axis
    perturbations.relativity = 1; 	% relativistic effects
    perturbations.srp = 1;          % solar radiation pressure
	perturbations.sun = 1;          % 3rd body gravity of Sun

    % ---------------
    % Initialization.
    % ---------------
    
    % initialize timing
    time = initialize_timing(UTC0_cal,t_duration,dt);
    data = initialize_data;
    
    % packages all structures into a single structure
    prop.data = data;
    prop.models = models;
    prop.perturbations = perturbations;
    prop.time = time;

    % TODO
	prop.epoch.GPS0_cal = time.GPS0_cal;
    prop.epoch.GPS0_MJD = time.GPS0_MJD;
    prop.epoch.GPS0_ws = time.GPS0_ws;
    prop.epoch.TT0_MJD = time.TT0_MJD;
    prop.epoch.t0 = time.t0;
    prop.epoch.UT10_MJD = time.UT10_MJD;
    prop.epoch.UTC0_cal = time.UTC0_cal;
    prop.epoch.UTC0_MJD = time.UTC0_MJD;

end