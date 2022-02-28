% mag_field  Calculates the magnetic field vector resolved in the ECI
% frame.
%
%   B_eci = mag_field(r,cotheta,phi,g,h,theta_gmst) returns the magnetic 
%   field vector "B_eci" resolved in the ECI frame given the geocentric 
%   distance "r", geocentric elevation "cotheta", geocentric longitude 
%   "phi", Gaussian-normalized spherical harmonic coefficients "g" and "h",
%   and the Greenwich mean sidereal time "theta_gmst".
%
% Copyright (c) 2021 Luke Neise
% Last Update: 2021-05-04



%% FUNCTION

% INPUT: r_eci - position vector in ECI frame [km]
%        cotheta - geocentric elevation (because theta = coelevation) [deg]
%                  NOTE: cotheta is essentially geocentric latitude
%        phi - geocentric longitude [deg]
%              NOTE: "phi" is used to match Wertz notation
%        g - Gaussian-normalized spherical harmonic coefficients [nT]
%        h - Gaussian-normalized spherical harmonic coefficients [nT]
%        theta_gmst - Greenwich mean sidereal time [rad]
% OUTPUT: B_eci - magnetic field vector resolved in the ECI frame [T]
function B_eci = mag_field(r_eci,cotheta,phi,g,h,theta_gmst)
    
    % radial distance
    r = norm(r_eci);
    
    % Earth mean equatorial radius [km]
    R_earth = 6378.1363;
    
    % coelevation [rad]
    theta = (90-cotheta);
    
    % converts coelevation and longtiude from deg to rad
    theta = theta*(pi/180);
    phi = phi*(pi/180);
    
    % compute Gauss function form of Legendre polynomials and its partial 
    % derivatives recursively
    P = zeros(5,5);
    dP = zeros(5,5);
    P(1,1) = 1;
    dP(1,1) = 0;
    for i = 2:5 %n=1,2,3,4
        for j = 1:5 %m=0,1,2,3,4
            n = i-1; m = j-1;
            if j > i
                continue;
            end
            if i == j
                P(i,j) = sin(theta)*P(i-1,j-1);
                dP(i,j) = sin(theta)*dP(i-1,j-1)+cos(theta)*P(i-1,j-1);
            elseif n == 1
                P(i,j) = cos(theta)*P(i-1,j);
                dP(i,j) = cos(theta)*dP(i-1,j)-sin(theta)*P(i-1,j);
            else
                K = ((n-1)^2-m^2)/(2*n-1)/(2*n-3);
                P(i,j) = cos(theta)*P(i-1,j)-K*P(i-2,j);
                dP(i,j) = cos(theta)*dP(i-1,j)-sin(theta)*P(i-1,j)-K*...
                    dP(i-2,j);
            end
        end
    end
    
    % gets rid of unnecessary n=0 dimension used for recursion
    P = P(2:end,:);
    dP = dP(2:end,:);
    
    % compute magnetic field strength in radial/south/east basis
    Br = 0;
    Bs = 0;
    Be = 0;
    for n = 1:4
        innerr = 0;
        inners = 0;
        innere = 0;
        for m = 0:n
            nidx = n;
            midx = m+1;
            innerr = innerr+(g(nidx,midx)*cos(m*phi)+h(nidx,midx)*sin(m*...
                phi))*P(nidx,midx);
            inners = inners+(g(nidx,midx)*cos(m*phi)+h(nidx,midx)*sin(m*...
                phi))*dP(nidx,midx);
            innere = innere+m*(-g(nidx,midx)*sin(m*phi)+h(nidx,midx)*...
                cos(m*phi))*P(nidx,midx);
        end
        Br = Br+((R_earth/r)^(n+2))*(n+1)*innerr;
        Bs = Bs+((R_earth/r)^(n+2))*inners;
        Be = Be+((R_earth/r)^(n+2))*innere;
    end
    Bs = Bs*(-1);
    Be = Be*(-1/sin(theta));
    
    % assigns cocoelivation to delta
    delta = cotheta*pi/180;
    
    % right ascension [rad]
    alpha = phi+theta_gmst;
    
    % ECI components of the magnetic field vector [nT]
    BI = (Br*cos(delta)+Bs*sin(delta))*cos(alpha)-Be*sin(alpha);
    BJ = (Br*cos(delta)+Bs*sin(delta))*sin(alpha)+Be*cos(alpha);
    BK = (Br*sin(delta)-Bs*cos(delta));
    
    % packages B vector and converts from nT to T
    B_eci = (1e-9)*[BI;BJ;BK];
    
end