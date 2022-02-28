%==========================================================================
%
% aero_torque  Aerodynamic torque acting on a satellite.
%
%   M_aero = aero_torque(r,v,A,rho,N_prin,b_prin,A_sat,CD)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   r       (3x1) [km] position vector resolved in ECI frame
%   v       (3x1) [km/s] inertial velocity vector resolved in ECI frame
%   A       (3x3) attitude matrix (ECI --> principal)
%   rho     (1x1) [kg/m^3] atmospheric density
%   N_prin  (Nx3) principal frame unit normal matrix
%   b_prin  (Nx3) [m] principal frame barycenter matrix
%   A_sat   (Nx1) [m^2] exposed surface area array
%   CD      (Nx1) drag coefficient array
%
% --------
% OUTPUTS:
% --------
%   M_aero	(3x1) [N.m] aerodynamic torque
%
% NOTE: "N_prin", "b_prin", "A_sat", and "CD" are ordered consistently and
% correspond to the surfaces of a component that is exposed to the
% environment. While a component may be exposed to the environment, one (or
% more) of its surfaces may not be; if this is the case, the corresponding
% row in "N_prin" is a 0 vector.
%
%==========================================================================
function M_aero = aero_torque(r,v,A,rho,N_prin,b_prin,A_sat,CD)

    % Earth rotational speed [rad/s]
    w_earth = 0.0000729211585530;
    
    % velocity of satellite w.r.t. incident flow [km/s]
    v_rel = v-cross([0;0;w_earth],r);
    
    % velocity of a satellite surface w.r.t. incident flow [km/s]
    V_eci = v_rel;
    
    % resolves V in principal frame and converts to m/s
    V = 1000*A*V_eci;
    
    % magnitude [m/s] of and unit vector in direction of V
    V_mag = norm(V);
    V_hat = V/V_mag;
    
    % initializes aerodynamic torque
    M_aero = [0;0;0];
    
    % calculates contribution of each surface to net aerodynamic torque
    for i = 1:length(A_sat)
        
        % extracts unit normal, center of pressure, area, and drag coeff.
        N = N_prin(i,:)';
        rs = b_prin(i,:)';
        Ai = A_sat(i);
        CDi = CD(i);

        % aerodynamic force acting on face [N]
        if idot(N,V_hat) > 0
            F_aero = -0.5*CDi*Ai*rho*V_mag^2*idot(V_hat,N)*V_hat;
        else
            F_aero = [0;0;0];
        end

        % adds contribution of surface to net aerodynamic torque [N.m]
        M_aero = M_aero+cross(rs,F_aero);
        
    end
    
end