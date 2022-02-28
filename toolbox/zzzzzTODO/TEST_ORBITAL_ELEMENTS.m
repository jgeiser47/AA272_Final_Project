clear;clc;close all;

% orbital elements
a = 10000;  % semi-major axis [km]
e = 0.4;    % eccentricity [-]
i = 40;     % inclination [deg]
Om = 25;    % right ascension of the ascending node [deg]
w = 80;     % argument of periapsis [deg]
nu0 = 60;    % true anomaly at epoch [deg]

% converts angles to radians
i = deg2rad(i);
Om = deg2rad(Om);
w = deg2rad(w);
nu0 = deg2rad(nu0);

% ECI position [km] and velocity [km/s] at epoch
[r_eci0,v_eci0] = oe2eci(a,e,i,Om,w,nu0);

% ECI position [km] throughout orbit
nu = (0:0.01:(2*pi))';
r_eci = zeros(3,length(nu));
for j = 1:length(nu)
    r_eci(:,j) = oe2eci(a,e,i,Om,w,nu(j));
end
rI = r_eci(1,:);
rJ = r_eci(2,:);
rK = r_eci(3,:);

% Earth gravitational parameter [km^3/s^2]
mu = 398600.4415;

% eccentricity vector at epoch
e_vec = ((norm(v_eci0)^2/mu)-(1/norm(r_eci0)))*r_eci0-(dot(r_eci0,...
    v_eci0)/mu)*v_eci0;

e_vec = a*e_vec;

% plot
figure;
hold on;
plot3(rI,rJ,rK,'linewidth',1.5);
plot3([0,e_vec(1)],[0,e_vec(2)],[0,e_vec(3)],'linewidth',1.5);
plot3([0,2*a],[0,0],[0,0],'color','k','linewidth',1.5);
plot3([0,0],[0,2*a],[0,0],'color','k','linewidth',1.5);
plot3([0,0],[0,0],[0,2*a],'color','k','linewidth',1.5);
hold off;
xlabel('$I\;[\mathrm{km}]$','interpreter','latex','fontsize',18);
ylabel('$J\;[\mathrm{km}]$','interpreter','latex','fontsize',18);
zlabel('$K\;[\mathrm{km}]$','interpreter','latex','fontsize',18);
grid on;
view(140,25);