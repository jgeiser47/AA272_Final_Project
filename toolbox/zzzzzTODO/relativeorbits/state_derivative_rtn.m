%==========================================================================
%
% state_derivative_rtn  Derivative of the relative RTN state of the deputy
% spacecraft with respect to the chief spacecraft.
%
%   dxdt = state_derivative_rtn(dx,d,ac)
%
% Author: Tamas Kis
% Last Update: 2021-08-16
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x_rel  	- (6×1 double) deputy relative state resolved in chief's RTN 
%             frame
%               --> 1-3. r_rel - relative position [m]
%               --> 4-6. v_rel - relative velocity [m/s]
%   d       - (3×1 double) differential perturbation acting on deputy,
%             resolved in chief's RTN frame [m/s^2]
%   u       - (3×1 double) control input [m/s^2]
%   ac      - (1×1 double) chief's semi-major axis [m]
%
% -------
% OUTPUT:
% -------
%   dx_rel  - (6×1 double) relative RTN state derivative
%               --> (1-3): relative velocity [m/s]
%               --> (4-6): relative acceleration [m/s^2]
%
%==========================================================================
function dx_rel = state_derivative_rtn(x_rel,d,u,ac)

    % unpacks relative RTN state [m][m/s^2]
    x = x_rel(1);
    y = x_rel(2);
    z = x_rel(3);
    xdot = x_rel(4);
    ydot = x_rel(5);
    zdot = x_rel(6);
    
    % unpacks differential perturbation [m/s^2]
    dx = d(1);
    dy = d(2);
    dz = d(3);
    
    % unpacks control input [m/s^2]
    ux = u(1);
    uy = u(2);
    uz = u(3);
    
 	% Earth gravitational parameter [m^3/s^2]
    mu_earth = 398600.4415e9;
    
    % chief's mean motion [rad/s]
    nc = sqrt(mu_earth/ac^3);
    
    % relative accelerations [m/s^2]
    xddot = 2*nc*ydot+nc^2*x-((mu_earth*(ac+x))/((ac+x)^2+y^2+z^2)^...
        (3/2))+(mu_earth/ac^2)+dx+ux;
    yddot = -2*nc*xdot-nc^2*y-((mu_earth*y)/((ac+x)^2+y^2+z^2)^(3/2))+...
        dy+uy;
    zddot = -((mu_earth*y)/((ac+x)^2+y^2+z^2)^(3/2))+dz+uz;
    
    % packages RTN state derivative [m/s][m/s^2]
    dx_rel = [xdot;ydot;zdot;xddot;yddot;zddot];

end