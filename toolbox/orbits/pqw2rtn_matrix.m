%==========================================================================
%
% pqw2rtn_matrix  Rotation matrix from perifocal (PQW) frame to 
% Hill/orbital (RTN) frame.
%
%   R_pqw2rtn = pqw2rtn_matrix(nu)
%
% Author: Tamas Kis
% Last Update: 2022-01-06
%
% REFERENCES:
%   [1] Vallado, "Fundamentals of Astrodynamics and Applications", 4th Ed.
%       (p. 169)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   nu          - (1×1 double) true anomaly [rad]
%
% -------
% OUTPUT:
% -------
%   R_pqw2rtn   - (3×3 double) rotation matrix (PQW --> RTN)
%
%==========================================================================
function R_pqw2rtn = pqw2rtn_matrix(nu)
    R_pqw2rtn = rot3(nu);
end