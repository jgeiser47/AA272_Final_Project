% TODO  TODO
%
% Author: Tamas Kis
% Last Update: 2021-08-02



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear;clc;close all;

% loads plot parameters and extracts needed parameters
pp = PLOT_PARAMETERS();



%% PARAMETERS

% epoch (UTC time) [YYYY,MM,DD,hh,mm,ss]
t0 = [2017,1,1,0,0,0];

% UT1/UTC offset [ms]
%   --> https://datacenter.iers.org/data/207/bulletinb-348.txt
dUT1 = 591.2975;

% orbital elements at epoch
a = 6853137;    % semi-major axis [m]
e = 0.001;      % eccentricity [-]
i = 0.9006;     % inclination [rad]
Om = 0;         % right ascension of the ascending node [rad]
w = 0;          % argument of periapsis [rad]
nu0 = pi/2;     % true anomaly [rad]

%Om = deg2rad(0);
%w = deg2rad(0);
%nu0 = deg2rad(0);

% Earth gravitational parameter [m^3/s^2]
mu = 398600.4415e9;

% mean motion [rad/s] and period [s]
n = a2n(a);
T = a2T(a);

% true latitude at epoch [rad]
u0 = w+nu0;

% time vector [s]
t = 0:1:2*T;



%% GEODETIC LATITUDE

% geodetic latitude from Kepler propagator
[~,~,~,~,~,~,~,lat,lon,~,theta_gmst] = kepler_propagator(t,a,e,i,Om,w,...
    nu0,t0,dUT1,mu);
lat_approx = rad2deg(asin(sin(u0+n*t)*sin(i)));

% error in the approximation
lat_error = lat_approx-lat;

% ----------------------------
% Plot of latitude with error.
% ----------------------------

% initializes figure
figure('position',pp.two_subplot_position);

% latitude plot
subplot(1,2,1);
hold on;
plot(t,lat_approx,'linewidth',pp.line_width,'color',pp.cardinal_red);
plot(t,lat,'k:','linewidth',pp.line_width);
hold off;
grid on;
xlabel('Time, $t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('Geodetic Latitude, $\phi\;[^{\circ}]$','interpreter',...
    'latex','fontsize',pp.axis_font_size);
legend('approximated latitude','true latitude','interpreter','latex',...
    'fontsize',pp.legend_font_size);

% latitude error plot
subplot(1,2,2);
plot(t,lat_error,'linewidth',pp.line_width,'color',pp.cardinal_red);
grid on;
xlabel('Time, $t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('Error in Latitude, $\Delta\phi\;[^{\circ}]$','interpreter',...
    'latex','fontsize',pp.axis_font_size);
legend("$("+scientific_notation_string(mean(lat_error))+")^{\circ}\pm"+...
    std(lat_error)+"^{\circ}$",'interpreter','latex','fontsize',...
    pp.legend_font_size);



%% GEODETIC LONGITUDE

% argument of latitude
u = u0+n*t;

% geodetic longitude
lon_approx = zeros(size(t));
Y = zeros(size(t));
for j = 1:length(t)
    numerator = cos(theta_gmst(j)-Om)*cos(u(j))+sin(theta_gmst(j)-Om)*sin(u(j))*cos(i);
    denominator = sqrt(1-sin(u(j))^2*sin(i)^2);
    Y(j) = -sin(theta_gmst(j)-Om)*cos(u(j))+cos(theta_gmst(j)-Om)*sin(u(j))*cos(i);
    if Y(j) > 0
        lon_approx(j) = rad2deg(acos(numerator/denominator));
    else
        lon_approx(j) = wrapTo180(360-rad2deg(acos(numerator/denominator)));
    end
end

% error in the approximation
lon_error = lon_approx-lon;

% deletes "jumps" from longitude error
[t_error,lon_error] = delete_data_points(t,lon_error,abs(lon_error)>1);

% -----------------------------
% Plot of longitude with error.
% -----------------------------

% initializes figure
figure('position',pp.two_subplot_position);

% longitude plot
subplot(1,2,1);
hold on;
plot(t,lon_approx,'linewidth',pp.line_width,'color',pp.cardinal_red);
plot(t,lon,'k:','linewidth',pp.line_width);
hold off;
grid on;
xlabel('Time, $t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('Geodetic Longitude, $\lambda\;[^{\circ}]$','interpreter',...
    'latex','fontsize',pp.axis_font_size);
legend('approximated longitude','true longitude','interpreter','latex',...
    'fontsize',pp.legend_font_size);

% longitude error plot
subplot(1,2,2);
plot(t_error,lon_error,'linewidth',pp.line_width,'color',pp.cardinal_red);
grid on;
xlabel('Time, $t\;[\mathrm{s}]$','interpreter','latex','fontsize',...
    pp.axis_font_size);
ylabel('Error in Longitude, $\Delta\lambda\;[^{\circ}]$','interpreter',...
    'latex','fontsize',pp.axis_font_size);
legend("$("+scientific_notation_string(mean(lon_error))+")^{\circ}\pm"+...
    std(lon_error)+"^{\circ}$",'interpreter','latex','fontsize',...
    pp.legend_font_size);