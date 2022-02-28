%==========================================================================
%
% gg_torque  Gravity gradient torque acting on a satellite.
%
%   M_gg = gg_torque(I,A,r)
%
% Author: Tamas Kis
% Last Update: 2021-09-09
%
% REFERENCES:
%   [1] D'Amico, "Gravity Gradient Torque, Stability, Damping", AA 279C 
%       Lecture 8 Slides
%   [2] Wertz, "Spacecraft Attitude Determination and Control" (pp. 566-
%       570)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   I       - (3×3 double) principal inertia tensor [kg.m^2]
%   A       - (3×3 double) attitude matrix (ECI --> principal) [-]
%   r       - (3×1 double) position resolved in ECI frame [m]
%
% -------
% OUTPUT:
% -------
%   M_gg    - (3×1 double) gravity gradient torque [N.m]
%
%==========================================================================
function M_gg = gg_torque(I,A,r)

    % Earth gravitational parameter [m^3/s^2]
    mu_earth = 398600.4415e9;
    
    % position vector magnitude [m]
    r_mag = norm(r);
    
    % radial basis vector
    R = r/r_mag;
    
    % c vector
    c = A*R;
    
    % extracts components of c vector
    cx = c(1);
    cy = c(2);
    cz = c(3);
    
    % moments of inertia about principal axes [kg.m^2]
    Ix = I(1,1);
    Iy = I(2,2);
    Iz = I(3,3);
    
    % gravity gradient torque [N.m]
    M_gg = (3*mu_earth/r_mag^3)*[(Iz-Iy)*cy*cz;
                                 (Ix-Iz)*cz*cx;
                                 (Iy-Ix)*cx*cy];

end