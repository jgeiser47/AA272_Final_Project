%==========================================================================
%
% eul2dcm313  3-1-3 Euler angles to direction cosine matrix (DCM).
%
%   A = eul2dcm313(phi,theta,psi)
%
% Authors: Tamas Kis, Luke Neise
% Last Update: 2021-09-09
%
% REFERENCES:
%   [1] Wertz, "Spacecraft Attitude Determination and Control", Eq. (E-14) 
%       (p. 763)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   phi     - (1×1 double) 3-1-3 Euler angle [rad]
%   theta   - (1×1 double) 3-1-3 Euler angle [rad]
%   psi     - (1×1 double) 3-1-3 Euler angle [rad]
%
% -------
% OUTPUT:
% -------
%   A       - (3×3 double) direction cosine matrix [-]
%
%==========================================================================
function A = eul2dcm313(phi,theta,psi)

    % precomputes trigonometric functions
    sp = sin(psi);
    cp = cos(psi);
    so = sin(theta);
    co = cos(theta);
    sf = sin(phi);
    cf = cos(phi);
     
    % defines direction cosine matrix
    A = [ cp*cf-sp*co*sf    cp*sf+sp*co*cf   sp*so;
         -sp*cf-cp*co*sf   -sp*sf+cp*co*cf   cp*so;
          so*sf            -so*cf            co];
 
end