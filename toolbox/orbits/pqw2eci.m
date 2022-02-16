%==========================================================================
%
% pqw2eci  PQW position and velocity to ECI position and velocity.
%
%   [r_eci,v_eci] = pqw2eci(r_pqw,v_pqw,i,Om,w)
%
% Author: Tamas Kis
% Last Update: 2022-01-06
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   r_pqw   - (3×1 double) position resolved in PQW frame
%   v_pqw   - (3×1 double) inertial velocity resolved in PQW frame
%   i       - (1×1 double) inclination [rad]
%   Om      - (1×1 double) right ascension of the ascending node [rad]
%   w       - (1×1 double) argument of periapsis [rad]
%
% -------
% OUTPUT:
% -------
%   r_eci   - (3×1 double) position resolved in ECI frame
%   v_eci   - (3×1 double) inertial velocity resolved in ECI frame
%
% -----
% NOTE:
% -----
%   --> The position and velocity vectors can be input in any units, but
%       they MUST be consistent. The output position and velocity vectors
%       will be in the same units.
%
%==========================================================================
function [r_eci,v_eci] = pqw2eci(r_pqw,v_pqw,i,Om,w)

    % rotation matrix from PQW frame to ECI frame
    R_pqw2eci = pqw2eci_matrix(i,Om,w);
    
    % position and inertial velocity resolved in ECI frame
    r_eci = R_pqw2eci*r_pqw;
    v_eci = R_pqw2eci*v_pqw;
    
end