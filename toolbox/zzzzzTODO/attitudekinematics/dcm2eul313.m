%==========================================================================
%
% dcm2eul313  Direction cosine matrix (DCM) to 3-1-3 Euler angles.
%
%   [phi,theta,psi] = dcm2eul313(A)
%
% Authors: Tamas Kis, Luke Neise
% Last Update: 2021-08-26
%
% REFERENCES:
%   [1] Wertz, "Spacecraft Attitude Determination and Control", Eq. 
%       (E-14) (p. 763)
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   A       - (3×3 double) direction cosine matrix
%
% --------
% OUTPUTS:
% -------- 
%   phi   	- (1×1 double) 3-1-3 Euler angle [rad]
%   theta   - (1×1 double) 3-1-3 Euler angle [rad]
%   psi     - (1×1 double) 3-1-3 Euler angle [rad]
%
% -----
% NOTE:
% -----
%   --> A33 is sometimes slightly greater than 1 or slightly less than 
%       negative 1, which would lead to an imaginary solution with the acos
%       function. Thus, we limit A33 to be in the domain [-1,1].
%   --> Alternatively, we could calculate theta directly without
%       restricting the domain, and then keep only the real portion of the
%       result.
%   --> Either of the above two methods works. However, we use the second
%       method as the first can still result in an error in Simulink.
%
%==========================================================================
function [phi,theta,psi] = dcm2eul313(A)
    phi = mod(atan2(A(3,1),-A(3,2)),2*pi);
    theta = acos(max(min(A(3,3),1),-1));
    psi = mod(atan2(A(1,3),A(2,3)),2*pi);    
end