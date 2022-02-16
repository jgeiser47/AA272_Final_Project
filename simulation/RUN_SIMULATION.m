%% RUN_SIMULATION.m
% Astrodynamics Toolbox
%
% TODO
%
% Author: Tamas Kis
% Last Update: 2022-02-07



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to all "Astrodynamics Toolbox" functions
addpath(genpath("../toolbox"));

% add path to external toolboxes
%addpath(genpath('../../MATLAB/toolboxes/State_Estimation_Toolbox/toolbox'));

% generate new orbits? (1 = yes, 0 = no)
generate_new_orbits = 1;



%% LOAD SETTINGS (src/settings)


% spacecraft settings
chief = set_chief;              % low drag attitude
deputy1 = set_deputy1(chief);   % high drag attitude
deputy2 = set_deputy2(chief);   % high drag attitude

models.density = 'Harris-Priester';
models.grav_mu = MU_EARTH;
models.grav_N = 0;
models.grav_R = R_EARTH;

perturb.drag = true;
perturb.relativity = false;
perturb.srp = true;

UTC_start = [2017,1,1,0,0,0];
duration = 12;
time_step = 10;
integrator = 'ABM8';
sim_start = 0;

prop = init_prop(models,perturb,UTC_start,duration,time_step,integrator,...
    sim_start);

chief_simdata = simulate_orbit(chief,prop);


%% RESULTS

% extracts chief orbit
X = chief_simdata.r_eci(1,:);
Y = chief_simdata.r_eci(2,:);
Z = chief_simdata.r_eci(3,:);

% extracts chief true orbit
%X_true = x_true(1,:);
%Y_true = x_true(2,:);
%Z_true = x_true(3,:);

figure;
hold on;
planet3D('earth cloudy');
plot3(X,Y,Z);
hold off;

figure;
hold on;
plot(chief_simdata.t,chief_simdata.g_eci(1,:));
plot(chief_simdata.t,chief_simdata.g_eci(2,:));
plot(chief_simdata.t,chief_simdata.g_eci(3,:));
hold off;

% prop = set_propagator;
% 
% 
% 
% %% SIMULATION DATA
% 
% % simulates new orbits
% if generate_new_orbits
%     
%     % simulates orbits
%     chief_simdata = simulate_orbit(chief,prop);
%     deputy1_simdata = simulate_orbit(deputy1,prop);
%     deputy2_simdata = simulate_orbit(deputy2,prop);
%     
%     % saves data into .mat files
%     save('simdata/chief_simdata.mat','chief_simdata');
%     save('simdata/deputy1_simdata.mat','deputy1_simdata');
%     save('simdata/deputy2_simdata.mat','deputy2_simdata');
% 
% % loads previously simulated orbits
% else
%     chief_simdata = struct2array(load('simdata/chief_simdata.mat'));
%     deputy1_simdata = struct2array(load('simdata/deputy1_simdata.mat'));
%     deputy2_simdata = struct2array(load('simdata/deputy2_simdata.mat'));
%     
% end
% 
% % combines all simulation data into one structure
% simdata.chief_simdata = chief_simdata;
% simdata.deputy1_simdata = deputy1_simdata;
% simdata.deputy2_simdata = deputy2_simdata;
% 
% 
% 
% %% STATE ESTIMATION
% 
% % simulates ground truth state
% x_true = true_state(prop,chief_simdata,deputy1_simdata,deputy2_simdata);
% 
% % initialize filter parameters
% filter = initialize_filter(x_true,chief,deputy1,deputy2,prop);
% 
% % simulates ground truth
% y = true_measurement(x_true,simdata,prop,filter);
% 
% % zero control input [m/s^2]
% u = zeros(3,prop.time.N);
% 
% % runs filter
% [x,P,tsol,rank_Ob,z_pre,z_post] = EKF(filter.f,filter.h,filter.F,...
%     filter.H,filter.Q,filter.R,prop.time.t,u,y,filter.x0,filter.P0,true);
% % [x,P,tsol,z_pre,z_post] = UKF(filter.f,filter.h,filter.Q,filter.R,...
% %    prop.time.t,u,y,filter.x0,filter.P0,true);
% 
% 
% 
% %% PROCESS NOISE COVARIANCE
% 
% Q = process_noise_covariance(x_true,filter.f,u);
% 
% 
% 
% %% RESULTS
% 
% % extracts chief orbit
% X = x(1,:);
% Y = x(2,:);
% Z = x(3,:);
% 
% % extracts chief true orbit
% X_true = x_true(1,:);
% Y_true = x_true(2,:);
% Z_true = x_true(3,:);
% 
% % figure;
% % hold on;
% % planet3D('earth cloudy');
% % plot3(X,Y,Z);
% % hold off;
% 
% 
% 
% %figure;plot(prop.time.t,X-X_true);title('x absolute error');
% %figure;plot(prop.time.t,Y-Y_true);title('y absolute error');
% %figure;plot(prop.time.t,Z-Z_true);title('z absolute error');
% 
% %figure;plot(prop.time.t,abs((X-X_true)./X_true));title('x relative error');
% %figure;plot(prop.time.t,abs((Y-Y_true)./Y_true));title('y relative error');
% %figure;plot(prop.time.t,abs((Z-Z_true)./Z_true));title('z relative error');
% 
% %figure;plot(prop.time.t,Y_true);title('y');
% 
% % extracts density results
% rhom1 = x(19,:);
% rhoM1 = x(20,:);
% rhom2 = x(21,:);
% rhoM2 = x(22,:);
% 
% % extracts true density results
% rhom1_true = x_true(19,:);
% rhoM1_true = x_true(20,:);
% rhom2_true = x_true(21,:);
% rhoM2_true = x_true(22,:);
% 
% h1 = 460;
% h2 = 500;
% 
% rho_est = zeros(size(rhoM1));
% % calculates chief density
% for i = 1:size(x,2)
%     UT1_MJD = prop.time.UT1_MJD(i);
%     r_eci = x(1:3,i);
%     v_eci = x(4:6,i);
%     GMST = UT1toGMST(UT1_MJD);
%     [r_ecef,v_ecef] = eci2ecef(r_eci,v_eci,GMST);
%     rho_est(i) = density_harris_priester_filter(r_eci,r_ecef,h1,...
%         rhom1(i),rhoM1(i),h2,rhom2(i),rhoM2(i),UT1_MJD);
% end
% 
% figure;
% hold on;
% plot(prop.time.t,chief_simdata.rho);
% plot(prop.time.t,rho_est);
% grid on;
% xlabel('time, s');
% ylabel('density, kg/m^3');
% legend('ground truth','estimate');
% hold off;
% 
% figure;
% hold on;
% plot(prop.time.t,rhom1);
% plot(prop.time.t,rhoM1);
% hold off;
% 
% figure;
% hold on;
% plot(prop.time.t,rhom2);
% plot(prop.time.t,rhoM2);
% hold off;
% 
% % calculates errors
% %rhom1_error = rhom1-rhom1_true;
% %rhoM1_error = rhoM1-rhoM1_true;
% %rhom2_error = rhom2-rhom2_true;
% %rhoM2_error = rhoM2-rhoM2_true;
% 
% 
% [lower_bound,upper_bound] = covariance_bounds(x,P);
% 
% rhom1_lower = lower_bound(19,:);
% rhoM1_lower = lower_bound(20,:);
% rhom2_lower = lower_bound(21,:);
% rhoM2_lower = lower_bound(22,:);
% 
% rhom1_upper = upper_bound(19,:);
% rhoM1_upper = upper_bound(20,:);
% rhom2_upper = upper_bound(21,:);
% rhoM2_upper = upper_bound(22,:);
% 
% opts.xunits = "g/km^{3}";
% opts.shaded = true;
% opts.dots = false;
% opts.color = [0,0.4470,0.7410];
% plot_filter_results(prop.time.t,rhom1,rhom1_true,rhom1_lower,rhom1_upper,opts)
% 
% % figure;
% % plot(prop.time.t,chief_simdata.h);
% % 
% % figure;
% % plot(prop.time.t,deputy1_simdata.h);
% % 
% % figure;
% % plot(prop.time.t,deputy2_simdata.h);
% % 
% % figure;
% % planet3D('earth cloudy');
% % hold on;
% % plot3(chief_simdata.r_eci(1,:),chief_simdata.r_eci(2,:),...
% %     chief_simdata.r_eci(3,:));
% % hold off;
% % 
% % figure;
% % ground_track(chief_simdata.lat,chief_simdata.lon)
% 
% 
% 
% 
% 
% 
% % 
% % 
% % 
% % 
% % 
% % %% PLOTS
% % 
% % close all;
% % 
% % figure;
% % hold on;
% % plot(t,x_true(1,:)-x(1,:));
% % hold off;
% % 
% % figure;
% % hold on;
% % plot(t,x_true(2,:)-x(2,:));
% % hold off;
% % 
% % figure;
% % hold on;
% % plot(t,x_true(3,:)-x(3,:));
% % hold off;
% % 
% % figure;
% % hold on;
% % plot(t,x_true(4,:)-x(4,:));
% % hold off;
% % 
% % figure;
% % hold on;
% % plot(t,x_true(5,:)-x(5,:));
% % hold off;
% % 
% % figure;
% % hold on;
% % plot(t,x_true(6,:)-x(6,:));
% % hold off;
% % 
% % figure;
% % hold on;
% % plot(t,x_true(7,:));
% % plot(t,x(7,:));
% % hold off;
% % 
% % % close plots if only running this section
% % %close all;
% % 
% % % produces plots
% % %make_plots(out);