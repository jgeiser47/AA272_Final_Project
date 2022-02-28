% euler_eqns_torque_free  Euler equations for torque-free motion.
%
%   dw = euler_eqns(t,w,I) returns the angular velocity derivative "dw" 
%   resolved in the principal frame given the angular velocity "w" (also 
%   resolved in the principal frame) and the principal inertia tensor "I".
%   Note that time, "t", is only included to make this function compatible
%   for use with ODE solvers.
%
% Copyright (c) 2021 Tamas Kis, Luke Neise
% Last Update: 2021-04-21



%% FUNCTION

% INPUT: t - time [s]
%        w - angular velocity vector resolved in the principal frame 
%            [rad/s]
%        I - principal inertia tensor
% OUTPUT: dw - angular velocity derivative resolved in the principal frame
%              [rad/s^2]
function dw = euler_eqns_torque_free(t,w,I)
    
    % angular velocities about principal axes
    wx = w(1);
    wy = w(2);
    wz = w(3);
    
    % moments of inertia about principal axes
    Ix = I(1,1);
    Iy = I(2,2);
    Iz = I(3,3);
    
    % angular velocity derivatives
    dwx = ((Iy-Iz)*wy*wz)/Ix;
    dwy = ((Iz-Ix)*wz*wx)/Iy;
    dwz = ((Ix-Iy)*wx*wy)/Iz;
    
    % packages angular velocity derivative vector
    dw = [dwx;dwy;dwz];

end