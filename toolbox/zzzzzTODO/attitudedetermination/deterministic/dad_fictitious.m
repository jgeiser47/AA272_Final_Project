%==========================================================================
%
% dad_fictitious  Deterministic attitude determination using fictitious 
% measurements.
%
%   A_body = dad_fictitious(M,U)
%
% Copyright (c) 2021 Tamas Kis, Luke Neise
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   M       (3x2) measurement matrix
%   U       (3x2) reference matrix
%
% --------
% OUTPUTS:
% --------
%   A_body  (3x3) body attitude matrix (reference --> measurement)
%
% -----
% NOTE:
% -----
%   --> Each column of "M" is a unit vector measurement resolved in the
%       measurement (body) frame
%   --> Each column of "U" is the corresponding ground truth reference unit
%       vector resolved in the reference (ECI) frame.
%
%==========================================================================
function A_body = dad_fictitious(M,U)

    % extracts measurement and reference vectors from M and U
    m1 = M(:,1);
    m2 = M(:,2);
    u1 = U(:,1);
    u2 = U(:,2);

    % fictitious measurements
    m1_fict = (m1+m2)/2;
    m1_fict = m1_fict/norm(m1_fict);
    m2_fict = (m1-m2)/2;
    m2_fict = m2_fict/norm(m2_fict);

    % fictitious reference vectors
    u1_fict = (u1+u2)/2;
    u1_fict = u1_fict/norm(u1_fict);
    u2_fict = (u1-u2)/2;
    u2_fict = u2_fict/norm(u2_fict);

    % assigns first fictitious measurement as reference axis
    pm = m1_fict;
    pu = u1_fict;

    % computes second measurement
    qm = cross(m1_fict,m2_fict)/norm(cross(m1_fict,m2_fict));
    qu = cross(u1_fict,u2_fict)/norm(cross(u1_fict,u2_fict));

    % completes the triad
    rm = cross(pm,qm);
    ru = cross(pu,qu);

    % assembles new M and U matrices
    M = [pm,qm,rm];
    U = [pu,qu,ru];

    % computes the DCM
    A_body = M*inv(U);

end