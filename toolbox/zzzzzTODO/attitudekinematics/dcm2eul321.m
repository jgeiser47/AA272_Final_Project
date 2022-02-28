%==========================================================================
%
% dcm2eul321  Direction cosine matrix (DCM) to 3-2-1 Euler angles.
%
%   [phi,theta,psi] = dcm2eul321(A)
%
% Author: Tamas Kis
% Last Update: 2021-08-26
%
% REFERENCES:
%   [1] Stevens, "Aircraft Control and Simulation", 3rd Ed., Eq. (1.3-11) 
%       (pp. 12)
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   A       - (3×3 double) direction cosine matrix [-]
%
% --------
% OUTPUTS:
% --------
%   phi 	- (1×1 double) 3-2-1 Euler angle [rad]
%   theta	- (1×1 double) 3-2-1 Euler angle [rad]
%   psi  	- (1×1 double) 3-2-1 Euler angle [rad]
%
% -----
% NOTE:
% -----
%   --> A13 is sometimes slightly greater than 1 or slightly less than 
%       negative 1, which would lead to an imaginary solution with the asin
%       function. Thus, we limit A13 to be in the domain [-1,1].
%   --> Alternatively, we could calculate theta directly without
%       restricting the domain, and then keep only the real portion of the
%       result.
%   --> Either of the above two methods works. However, we use the second
%       method as the first can still result in an error in Simulink.
%
%==========================================================================
function [phi,theta,psi] = dcm2eul321(A)
    phi = atan2(A(2,3),A(3,3));
    theta = -asin(max(min(A(1,3),1),-1));
    psi = atan2(A(1,2),A(1,1));
end