%% FILTER.m TODO
% Astrodynamics Toolbox
%
% TODO
%
% Author: Tamas Kis
% Last Update: 2022-02-21



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to all "Astrodynamics Toolbox" functions
addpath(genpath('../toolbox'));

% adds path to all "State Space Discretization and Linearization Toolbox" 
% functions
addpath(genpath("../../State_Space_Discretization_and_Linearization_To"+...
    "olbox-MATLAB/toolbox"));

% adds path to all "State Estimation Toolbox" functions
addpath(genpath('../../State_Estimation_Toolbox/toolbox'));



%% FILTER SETTINGS

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
duration = 1.5;

% ------------------------------
% Process and measurement noise.
% ------------------------------

% position [m] and velocity [m/s] standard deviations
r_std = 5;
v_std = 0.005;

% measurement noise covariance
R = diag([r_std,r_std,r_std,v_std,v_std,v_std].^2);

% process noise covariance
%R = (1/100)^2*Q;
%Q = R/100;
Q = R;

% ------------------------------
% Parameters for discretization.
% ------------------------------

% time step [s]
dt = 10;

% integrator ('TODO')
integrator = 'RK4';



%% FILTERING

% initialize filter parameters
filter = initialize_filter(UTC_start,duration);

% ------------------------------
% Continuous-time system.
% ------------------------------

% continuous-time, nonlinear dynamics equation
f = @(x,u,t) dynamics_model(t,x,filter,sat);

% continuous-time, nonlinear measurement equation
h = @(x,t) measurement_model(t,x,filter,sat);


% load simulation data
simdata = load_simdata('swarm_ex_1/simdata/chief_simdata.mat');

% generate true states and measurements
x_true = true_state(simdata);
y = true_measurement(x_true,R,simdata,filter);

% time vector
t = simdata.t;

% discretize system
fd = f2fd_num(f,dt,[],integrator);
hd = h2hd_num(h,dt);

% Jacobians
F = @(x,u,k) error_stm_num(f,x,u,k2t(k,dt),dt);
H = @(x,k) hd2H_num(hd,x,k);

% initial conditions
x0 = sat.ECI;
P0 = Q*10;

[x,P,tsol,rank_Ob,z_pre,z_post] = EKF(fd,hd,F,H,Q,R,t,[],y,x0,P0,true);

%%



figure;
hold on;
plot(t(5:end),z_pre(1,5:end),'k.');
plot(t(5:end),z_post(1,5:end),'r.');
hold off;
legend('pre-fit residuals','post-fit residuals','interpreter','latex',...
    'fontsize',14,'location','best');
xlabel('time $[\mathrm{s}]$','interpreter','latex','fontsize',18);
ylabel('$r_{X}$ measurement residuals $[\mathrm{m}]$','interpreter',...
    'latex','fontsize',18);
grid on;

%% 

figure;
hold on;
plot(t(5:end),z_pre(2,5:end),'k.');
plot(t(5:end),z_post(2,5:end),'r.');
hold off;

figure;
hold on;
plot(t(5:end),z_pre(3,5:end),'k.');
plot(t(5:end),z_post(3,5:end),'r.');
hold off;



figure;
hold on;
plot(t(5:end),z_pre(4,5:end),'k.');
plot(t(5:end),z_post(4,5:end),'r.');
hold off;

figure;
hold on;
plot(t(5:end),z_pre(5,5:end),'k.');
plot(t(5:end),z_post(5,5:end),'r.');
hold off;

figure;
hold on;
plot(t(5:end),z_pre(6,5:end),'k.');
plot(t(5:end),z_post(6,5:end),'r.');
hold off;

%%
M = 2;


[lower_bound,upper_bound] = covariance_bounds(x,P,M);

x_lower = lower_bound(1,:);
y_lower = lower_bound(2,:);
z_lower = lower_bound(3,:);
vx_lower = lower_bound(4,:);
vy_lower = lower_bound(5,:);
vz_lower = lower_bound(6,:);

x_upper = upper_bound(1,:);
y_upper = upper_bound(2,:);
z_upper = upper_bound(3,:);
vx_upper = upper_bound(4,:);
vy_upper = upper_bound(5,:);
vz_upper = upper_bound(6,:);


opts.xunits = "m";
opts.shaded = true;
opts.dots = false;
opts.color = [0,0.4470,0.7410];
opts.error = true;
opts.M = M;

opts.name = '$r_{I}$';
opts.xunits = "m";
plot_filter_results(t',x(1,:),x_lower,x_upper,x_true(1,:),opts)

opts.name = '$r_{J}$';
opts.xunits = "m";
plot_filter_results(t',x(2,:),y_lower,y_upper,x_true(2,:),opts)

opts.name = '$r_{K}$';
opts.xunits = "m";
plot_filter_results(t',x(3,:),z_lower,z_upper,x_true(3,:),opts)

opts.name = '$v_{I}$';
opts.xunits = 'm/s';
plot_filter_results(t',x(4,:),vx_lower,vx_upper,x_true(4,:),opts)

opts.name = '$v_{J}$';
opts.xunits = 'm/s';
plot_filter_results(t',x(5,:),vy_lower,vy_upper,x_true(5,:),opts)

opts.name = '$v_{K}$';
opts.xunits = 'm/s';
plot_filter_results(t',x(6,:),vz_lower,vz_upper,x_true(6,:),opts)

%figure;
%plot(t,x(1,:)-x_true(1,:));