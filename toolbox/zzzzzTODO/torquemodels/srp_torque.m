%==========================================================================
%
% srp_torque  Solar radiation pressure (SRP) torque acting on a satellite.
%
%   M_srp = srp_torque(S_hat,eclipse,N_prin,b_prin,A_sat,Cd,Cs)
%
% Author: Tamas Kis
% Last Update: 2021-09-11
%
% REFERENCES:
%   [1] D'Amico, "Other Perturbation Torques", AA 279C Lecture 9 Slides
%       (pp. 7-9)
%   [2] Wertz, "Spacecraft Attitude Determination and Control" (pp. 570-
%       573)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   S_hat   - (3×1 double) direction of Sun (unit vector pointing from 
%             satellite towards Sun) resolved in principal frame
%   eclipse - (logical) true if satellite is in eclipse, false if not
%   N_prin	- (N×3 double) principal frame unit normal matrix
%   b_prin	- (N×3 double) principal frame barycenter matrix [m]
%   A_sat   - (N×1 double) exposed surface area array [m^2]
%   Cd   	- (N×1 double) coefficient of diffuse reflection array
%   Cs    	- (N×1 double) coefficient of specular reflection array
%
% -------
% OUTPUT:
% -------
%   M_srp   - (3×1 double) solar radiation pressure torque [N.m]
%
% -----
% NOTE:
% -----
%   --> "N_prin", "b_prin", "A_sat", "Cd", and "Cs" are ordered 
%       consistently and correspond to the surfaces of a component that is 
%       exposed to the environment. While a component may be exposed to the 
%       environment, one (or more) of its surfaces may not be; if this is 
%       the case, the corresponding row in "N_prin" is a 0 vector.
%
%==========================================================================
function M_srp = srp_torque(S_hat,eclipse,N_prin,b_prin,A_sat,Cd,Cs)
    
    % speed of light [m/s]
    c = 2.9979e8;
    
    % mean integrated solar energy flux [W/m^2]
    Fe = 1358;
    
    % solar radiation pressure [N/m^2]
    P = Fe/c;
    
    % initializes aerodynamic torque
    M_srp = [0;0;0];
    
    % calculates contribution of each surface to net SRP torque
    for i = 1:length(A_sat)

        % extracts unit normal, center of pressure, area, and surface
        % coefficients
        N = N_prin(i,:)';
        rs = b_prin(i,:)';
        Ai = A_sat(i);
        Cdi = Cd(i);
        Csi = Cs(i);

        % solar radiation pressure force acting on face [N]
        if (idot(S_hat,N) > 0) && (~eclipse)
            F_srp = -P*Ai*((1-Csi)*S_hat+2*(Csi*idot(S_hat,N)+Cdi/3)*...
                N)*idot(S_hat,N);
        else
            F_srp = [0;0;0];
        end

        % adds contribution of surface to net SRP torque [N.m]
        M_srp = M_srp+cross(rs,F_srp);

    end
    
end