%==========================================================================
%
% quat2dcm  Quaternion to direction cosine matrix (DCM).
%
%   A = quat2dcm(q)
%
% Authors: Tamas Kis, Luke Neise
% Last Update: 2021-08-26
%
% REFERENCES:
%   [1] Wertz, "Spacecraft Attitude Determination and Control", Eq. 
%       (12-13a) (p. 414)
%   [2] Wertz, "Spacecraft Attitude Determination and Control", Eq. (E-8) 
%       (p. 762)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   q   - (4×1 double) quaternion [-]
%
% -------
% OUTPUT:
% -------
%   A   - (3×3 double) direction cosine matrix [-]
%
%==========================================================================
function A = quat2dcm(q)

    % unpacks quaternion
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);

    % direction cosine matrix
    A = [q1^2-q2^2-q3^2+q4^2    2*(q1*q2+q3*q4)       2*(q1*q3-q2*q4);
         2*(q1*q2-q3*q4)       -q1^2+q2^2-q3^2+q4^2   2*(q2*q3+q1*q4);
         2*(q1*q3+q2*q4)        2*(q2*q3-q1*q4)      -q1^2-q2^2+q3^2+q4^2];

end