%==========================================================================
%
% gyro_actuator_command  Angular accelerations (i.e. commands to gyroscopic
% actuators) required to produce desired control torque.
%
%   dw = gyro_actuator_command(w,Mc,Aw,Iw)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   w       (nwx1) [rad/s] angular velocity of each gyroscopic actuator
%   Mc      (3x1) [N.m] control torque
%   Aw      (3xnw) mounting matrix
%   Iw      (1x1) [kg.m^2] moment of inertia of each gyroscopic actuator
%
% --------
% OUTPUTS:
% --------
%   dw      (nwx1) [rad/s^2] angular acceleration required of each 
%                            gyroscopic actuator to produce Mc
%
%==========================================================================
function dw = gyro_actuator_command(w,Mc,Aw,Iw)
   
    % psuedoinverse of Aw
    Aw_star = Aw'*inv(Aw*Aw');
    
    % angular momenta of reaction wheels [kg.m^2/s]
    Lw = Iw*w;
    
    % required angular acceleration [rad/s^2]
    dw = Aw_star*(-Mc)*(1/Iw);

end