%==========================================================================
%
% drag  Perturbing acceleration due to atmospheric drag resolved in the ECI
% frame.
%
%   f_D_eci = drag(v_ecef,R_ecef2eci,rho,B)
%
% Author: Tamas Kis
% Last Update: 2022-02-01
%
% REFERENCES:
%   [1] Montenbruck and Gill, "Satellite Orbits" (pp. 83-86)
%   [2] Vallado, "Fundamentals of Astrodynamics and Applications", 4th Ed.
%       (pp. 551-552)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   v_ecef      - (3×1 double) ECEF velocity resolved in ECEF frame [m/s]
%   R_ecef2eci  - (3×3 double) rotation matrix (ECEF --> ECI)
%   rho         - (1×1 double) atmospheric density [kg/m^3]
%   B           - (1×1 double) ballistic coefficient [m^2/kg]
%
% -------
% OUTPUT:
% -------
%   f_D_eci     - (3×1 double) perturbing acceleration due to atmospheric 
%                 drag resolved in ECI frame [m/s^2]
%
%==========================================================================
function f_D_eci = drag(v_ecef,R_ecef2eci,rho,B)

    % relative velocity of satellite with respect to the atmosphere, 
    % resolved in the ECEF frame [m/s]
    v_rel_ecef = v_ecef;

    % relative velocity of satellite with respect to the atmosphere, 
    % resolved in the ECI frame [m/s]
    v_rel = R_ecef2eci*v_rel_ecef;
    
    % acceleration due to atmospheric drag [m/s]
    f_D_eci = -0.5*B*rho*inorm(v_rel)*v_rel;

end