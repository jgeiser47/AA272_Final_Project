%==========================================================================
%
% pqw2rtn  PQW position and velocity to RTN position and velocity.
%
%   [r_rtn,v_rtn] = pqw2rtn(r_pqw,v_pqw,nu)
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-06-01
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   r_pqw   (3x1) [km] position vector resolved in perifocal frame
%   v_pqw   (3x1) [km/s] inertial velocity vector resolved in perifocal 
%                        frame
%   nu      (1x1) [rad] true anomaly
%
% --------
% OUTPUTS:
% --------
%   r_rtn  	(3x1) [km] position vector resolved in RTN frame
%   v_rtn   (3x1) [km/s] inertial velocity vector resolved in RTN frame
%
%==========================================================================
function [r_rtn,v_rtn] = pqw2rtn(r_pqw,v_pqw,nu)
    
    % defines rotation matrix
    R_pqw2rtn = [ cos(nu)   sin(nu)   0;
                 -sin(nu)   cos(nu)   0;
                  0         0         1];
    
    % position [km] and inertial velocity [km/s] resolved in RTN frame
    r_rtn = R_pqw2rtn*r_pqw;
    v_rtn = R_pqw2rtn*v_pqw;
    
end