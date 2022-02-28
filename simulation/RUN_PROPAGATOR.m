%% RUN_PROPAGATOR.m
% Astrodynamics Toolbox
%
% Runs an orbit propagator to simulate an orbit.
%
% Author: Tamas Kis
% Last Update: 2022-02-21



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to all "Astrodynamics Toolbox" functions
addpath(genpath('../toolbox'));



%% PROPAGATOR SETTINGS

% -----------
% File paths.
% -----------

% path to folder storing settings
simulation_path = 'swarm_ex_1';
addpath(genpath(simulation_path));

% file name for simulation data
file_name = 'chief_simdata';

% --------------------
% Satellite selection.
% --------------------

sat = chief_parameters;

% -----
% Time.
% -----

% initial UTC [y,mo,d,h,m,s]
UTC_start = [2017,1,1,0,0,0];

% simulation duration [h]
%duration = 1.5;
duration = 6;

% --------------------------------
% Integator (ODE solver) settings.
% --------------------------------

% time step [s]
time_step = 10;

% integrator ('RK4' or 'ABM8')
%integrator = 'ABM8';
integrator = 'RK4';

% simulation start time [s]
sim_start = 0;

% -------
% Models.
% -------

% density model ('Exponential', 'Harris-Priester', 'Jacchia-Bowman 2008',
% 'Jacchia-Roberts', 'NRLMSISE-00', or 'NRLMSISE-00 MATLAB')
models.density = 'NRLMSISE-00';

% maximum degree/order for gravity model
models.grav_N = 120;

% --------------
% Perturbations.
% --------------

perturb.drag = true;
perturb.relativity = true;
perturb.srp = true;
perturb.moon = true;
perturb.sun = true;



%% SIMULATION

% initialize propagator
prop = initialize_propagator(models,perturb,UTC_start,duration,...
    time_step,integrator,sim_start);

% run simulation
simdata = simulate_orbit(sat,prop);



%% SAVING DATA

% creates "simdata" directory if needed
simdata_file_path = strcat(simulation_path,'/simdata/');
if ~exist(simdata_file_path,'dir')
    mkdir(simdata_file_path)
end

% saves simulation data
save(strcat(simdata_file_path,file_name,'.mat'),'simdata');