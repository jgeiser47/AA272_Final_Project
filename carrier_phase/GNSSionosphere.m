% s3_GNSSionosphere.m
%
% -------------------------------------------------------------------------
% USAGE: I = s3_GNSSionosphere
% -------------------------------------------------------------------------
% DESCRIPTION: Computes the Ionosphere path delay using the Klobuchar model
% -------------------------------------------------------------------------
% INPUTS:  GNSSt  - GNSS time of the observation 
%          r_rcv  - Position (or state) of receiver at GNSSt [m]
%          r_sat  - Position of (or state) the GNSS satellite at time of
%                   transmission [m]
%          Kalpha - Klobuchar alpha parameters
%          Kbeta  - Klobuchar beta parameters
% -------------------------------------------------------------------------
% OUTPUTS: I - The GNSS ionospheric path delay in [m]
% -------------------------------------------------------------------------
% REFERENCE: "Implementation of GNSS Ionospheric Models in gLAB"
%            by Deimos Ibanez Segura
% -------------------------------------------------------------------------
% SEE ALSO:
% -------------------------------------------------------------------------
% AUTHOR: deddy
%         Jul 21, 2016
% -------------------------------------------------------------------------
% COPYRIGHT:
% Built using S3: SLAB Satellite Software, (c) 2016
% https://slab.stanford.edu/s3/
% -------------------------------------------------------------------------

function I = GNSSionosphere(GPSsecond,r_rcv,r_sat,Kalpha,Kbeta)

r_rcv = r_rcv(:);
r_sat = r_sat(:);

r_rcv = r_rcv(1:3);
r_sat = r_sat(1:3);

% Declare Constants
R_EARTH = 6378137.0;   % [m] WGS84 value
C_LIGHT = 299792458.0; % [m/s]

% Azimuth and Elevation of GNSS sat with respect to user
rho_enu = ecef2enu(r_sat,r_rcv);
[A,E,~] = enu2aer(rho_enu);
A = deg2rad(A); E = deg2rad(E);

[LONu,LATu,~] = ecef2geod(r_rcv);
LONu = deg2rad(LONu); LATu = deg2rad(LATu);

% Convert To Semi-Circles
LONu = LONu/pi;
LATu = LATu/pi;

E = E/pi;

hiono = 350e3; % GPS Ionosphere height

%% BDS ICD Formulation

% A1 = 5e-9;
% A3 = 50400;
% 
% psi = pi/2 - E - asin(R_EARTH/(R_EARTH + hiono)*cos(E));
% 
% LATm = asin(sin(LATu)*cos(psi) + cos(LATu)*sin(psi)*cos(A));
% LONm = LONu + asin(sin(psi)*sin(A)/cos(LATm));
% 
% t = mod(43200*LONm + GNSSt.second,86400);
% 
% if t > 86400
%     t = t - 86400;
% elseif t < 0
%     t = t + 86400;
% end
% 
% A2 = Kalpha(1)*abs(LATm)^0 + Kalpha(2)*abs(LATm)^1 + ...
%      Kalpha(3)*abs(LATm)^2 + Kalpha(4)*abs(LATm)^3;
% 
% if A2 < 0
%     A2 = 0;
% end
% 
% A4 = Kbeta(1)*abs(LATm)^0 + Kbeta(2)*abs(LATm)^1 + ...
%      Kbeta(3)*abs(LATm)^2 + Kbeta(4)*abs(LATm)^3;
%  
% if A4 >= 172800
%     A4 = 172800;
% elseif A4 < 72000
%     A4 = 72000;
% end
% 
% if abs(t-A3) < A4/4;
%     I = A1 + A2*cos(abs(2*pi*(t-A3)/A4));
% else
%     I = A1;
% end
% 
% I = I*C_LIGHT;

%% GPS Ionosphere Model

A1 = 5e-9;
A3 = 50400;

PSI = 0.0137/(E+0.11) - 0.022;

PHIl = LATu + PSI*cos(A);

if PHIl > 0.416
    PHIl = 0.416;
elseif PHIl < -0.416
    PHIl = -0.416;
end

LAMl = LONu + (PSI*sin(A))/cos(PHIl*pi);

PHIm = PHIl + 0.064*cos((LAMl-1.617)*pi);

t = mod(4.32e4*LAMl+GPSsecond,86400);

if t > 86400
    t = t-86400;
elseif t < 0
    t = t + 86400;
end

F = 1 + 16*(0.53-E)^3;

A2 = Kalpha(1)*(PHIm)^0 + Kalpha(2)*(PHIm)^1 + ...
     Kalpha(3)*(PHIm)^2 + Kalpha(4)*(PHIm)^3;
 
if A2 < 0
    A2 = 0;
end
 
A4 = Kbeta(1)*(PHIm)^0 + Kbeta(2)*(PHIm)^1 + ...
     Kbeta(3)*(PHIm)^2 + Kbeta(4)*(PHIm)^3;

if A4 >= 172800
    A4 = 172800;
elseif A4 < 72000
    A4 = 72000;
end

% Ionospheric Time Delay
x = 2*pi*(t-A3)/A4;

if abs(x) > 1.57
    T = F*A1;
else
    T = F*(A1 + A2*(1-x^2/2 + x^4/24));
end

I = T*C_LIGHT;

% %% Full Model
% 
% LONu = LONu*pi;
% LATu = LATu*pi;
% E    = E*pi;
% 
% A1 = 5e-9;
% A3 = 50400;
% 
% PSI = pi/2 - E - asin(R_EARTH/(R_EARTH + 350e3)*cos(E));
% 
% PHIl = asin(sin(LATu)*cos(PSI) + cos(LATu)*sin(PSI)*cos(A));
% LAMl = LONu + (PSI*sin(A))/cos(PHIl);
% 
% PHIm = asin(sin(PHIl)*sind(78.3) + cos(PHIl)*cosd(78.3)*cos(LAMl - 291.0*pi/180));
% 
% t = mod(4.32e4*LAMl+GNSSt.second,86400);
% 
% if t > 86400
%     t = t-86400;
% elseif t < 0
%     t = t + 86400;
% end
% 
% F = sqrt(1 - (R_EARTH/(R_EARTH+350e3)*cos(E))^2);
% 
% A2 = Kalpha(1)*(PHIm)^0 + Kalpha(2)*(PHIm)^1 + ...
%      Kalpha(3)*(PHIm)^2 + Kalpha(4)*(PHIm)^3;
%  
% if A2 < 0;
%     A2 = 0;
% end
%  
% A4 = Kbeta(1)*(PHIm)^0 + Kbeta(2)*(PHIm)^1 + ...
%      Kbeta(3)*(PHIm)^2 + Kbeta(4)*(PHIm)^3;
% 
% if A4 >= 172800
%     A4 = 172800;
% elseif A4 < 72000
%     A4 = 72000;
% end
% 
% % Ionospheric Time Delay
% x = 2*pi*(t-A3)/A4;
% 
% if abs(x) > pi/2
%     T = F*A1;
% else
%     T = F*(A1 + A2*(1-x^2/2 + x^4/24));
% end
% 
% I = T*C_LIGHT;

end
