%==========================================================================
%
% nominal_attitude  Attitude matrix (DCM) relating the principal frame to
% the ECI frame for the nominal spacecraft attitude.
%
%   A_nom = nominal_attitude(r,v,R_body2prin)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   r               (3x1) [km] position vector resolved in ECI frame
%   v               (3x1) [km/s] inertial velocity vector resolved in ECI 
%                                frame
%   A               (3x3) attitude matrix (ECI --> principal)
%   R_body2prin     (3x3) rotation matrix (body --> principal)
%
% --------
% OUTPUTS:
% --------
%   A_nom	(3x1) nominal attitude matrix (ECI --> principal)
%
%==========================================================================
function A_nom = nominal_attitude(r,v,R_body2prin)

    % obtains rotation matrix from ECI frame to RTN frame
    R_eci2rtn = eci2rtn_matrix(r,v);
    
    % rotation matrix from RTN frame to TNR frame
    R_rtn2tnr = [0   1   0;
                 0   0   1;
                 1   0   0];
    
    % nominal attitude matrix
    A_nom = R_body2prin*R_rtn2tnr*R_eci2rtn;
    
end