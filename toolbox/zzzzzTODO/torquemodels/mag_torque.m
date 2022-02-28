%==========================================================================
%
% mag_torque  Magnetic torque acting on a satellite.
%
%   M_mag = mag_torque(B_eci,A,R_body2prin,m_mag)
%
% Authors: Tamas Kis, Luke Neise
% Last Update: 2021-09-09
%
% REFERENCES:
%   [1] D'Amico, "Other Perturbation Torques", AA 279C Lecture 9 Slides
%       (p. 4)
%   [2] Wertz, "Spacecraft Attitude Determination and Control", Eq. 17-58 
%       (p. 575)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   B_eci   	- (3×1 double) magnetic field vector resolved in ECI frame 
%                 [T]
%   A           - (3×3 double) attitude matrix (ECI --> principal) [-]
%   R_body2prin - (3×3 double) rotation matrix (body --> principal) [-]
%   m_mag    	- (3×1 double) satellite magnetic dipole moment [A.m^2]
%
% -------
% OUTPUT:
% -------
%   M_mag       - (3×1 double) magnetic torque [N.m]
%
%==========================================================================
function M_mag = mag_torque(B_eci,A,R_body2prin,m_mag)

    % magnetic field vector resolved in principal frame [T]
    B = A*B_eci;

    % magnetic dipole moment resolved in body frame [A.m^2]
    m_body = [0;0;m_mag];
    
    % magnetic dipole moment resolved in principal frame [A.m^2]
    m = R_body2prin*m_body;

    % magnetic torque [N.m]
    M_mag = cross(m,B);

end