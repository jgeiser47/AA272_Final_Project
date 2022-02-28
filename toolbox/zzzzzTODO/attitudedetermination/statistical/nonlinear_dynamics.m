%==========================================================================
%
% nonlinear_dynamics  Discrete-time nonlinear dynamics equation.
%
% Propagates the state vector to the next time step using the discrete-time
% nonlinear dynamics equation.
%
%   xk1 = nonlinear_dynamics(xk,uk,dt,I)
%
% Copyright (c) 2021 Luke Neise, Tamas Kis
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   xk    	(7x1) state vector at current time step
%   uk      (3x1) [N.m] control input at current time step
%   dt     	(1x1) [s] time step
%   I   	(3x3) [kg.m^2] principal inertia tensor
%
% --------
% OUTPUTS:
% --------
%   xk1     (7x1) state vector at next time step
%
%==========================================================================
function xk1 = nonlinear_dynamics(xk,uk,dt,I)

    % unpacks state vector
    q1 = xk(1); 
    q2 = xk(2);
    q3 = xk(3);
    q4 = xk(4);
    wx = xk(5);
    wy = xk(6);
    wz = xk(7);
    
    % unpacks control input
    Mcx = uk(1);
    Mcy = uk(2);
    Mcz = uk(3);
    
    % unpacks inertia tensor
    Ix = I(1,1);
    Iy = I(2,2);
    Iz = I(3,3);

    % propagates state variables one time step 
    q1_new = q1+0.5*dt*(wz*q2-wy*q3+wx*q4);
    q2_new = q2+0.5*dt*(-wz*q1+wx*q3+wy*q4);
    q3_new = q3+0.5*dt*(wy*q1-wx*q2+wz*q4);
    q4_new = q4+0.5*dt*(-wx*q1-wy*q2-wz*q3);
    wx_new = wx+((Mcx+(Iy-Iz)*wy*wz)/Ix)*dt;
    wy_new = wy+((Mcy+(Iz-Ix)*wx*wz)/Iy)*dt;
    wz_new = wz+((Mcz+(Ix-Iy)*wx*wy)/Iz)*dt;

    % packages propagated state vector
    xk1 = [q1_new;q2_new;q3_new;q4_new;wx_new;wy_new;wz_new];

end