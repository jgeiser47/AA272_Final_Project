%==========================================================================
%
% eci2pqw  ECI position and velocity to PQW position and velocity.
%
%   [r_pqw,v_pqw] = eci2pqw(r_eci,v_eci,i,Om,w)
%
% Author: Tamas Kis
% Last Update: 2022-01-06
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   r_eci   - (3×1 double) position resolved in ECI frame
%   v_eci   - (3×1 double) inertial velocity resolved in ECI frame
%   i       - (1×1 double) inclination [rad]
%   Om      - (1×1 double) right ascension of the ascending node [rad]
%   w       - (1×1 double) argument of periapsis [rad]
%
% --------
% OUTPUTS:
% --------
%   r_pqw   - (3×1 double) position resolved in PQW frame
%   v_pqw   - (3×1 double) inertial velocity resolved in PQW frame
%
% -----
% NOTE:
% -----
%   --> The position and velocity vectors can be input in any units, but
%       they MUST be consistent. The output position and velocity vectors
%       will be in the same units.
%
%==========================================================================
function [r_pqw,v_pqw] = eci2pqw(r_eci,v_eci,i,Om,w)

    % rotation matrix from ECI frame to PQW frame
    R_eci2pqw = eci2pqw_matrix(i,Om,w);
    
    % position and inertial velocity resolved in PQW frame
    r_pqw = R_eci2pqw*r_eci;
    v_pqw = R_eci2pqw*v_eci;
    
end