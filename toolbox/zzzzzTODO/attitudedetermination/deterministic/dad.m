%==========================================================================
%
% dad  Deterministic attitude determination.
%
%   A = dad(M,U)
%
% Copyright (c) 2021 Tamas Kis, Luke Neise
% Last Update: 2021-05-31
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   M       (3xN) measurement matrix
%   U       (3xN) reference matrix
%
% --------
% OUTPUTS:
% --------
%   A_body   (3x3) body attitude matrix (ECI --> body)
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
function A = dad(M,U)

    % determines number of measurements
    N = size(M,2);
    
    % DCM for 1 measurement
    if N == 1
        
        % measurement and reference vectors
        m = M(1:3,1);
        u = U(1:3,1);
        
        % computes DCM
        G = [idot(u,m)          -norm(cross(u,m))   0;
             norm(cross(u,m))    idot(u,m)          0;
             0                   0                  1];
        F = [u    (m-idot(u,m)*u)/norm(m-idot(u,m)*u)    cross(m,u)];
        A = F*G*inv(F);

    % DCM for 2 measurements
    elseif N == 2
        
        % extracts measurement and reference vectors from M and U
        m1 = M(:,1);
        m2 = M(:,2);
        u1 = U(:,1);
        u2 = U(:,2);
        
        % assigns the dominant measurement to the reference axis
        pm = m1;
        pu = u1;

        % computes second measurement
        qm = cross(m1,m2)/norm(cross(m1,m2));
        qu = cross(u1,u2)/norm(cross(u1,u2));

        % completes the triad
        rm = cross(pm,qm);
        ru = cross(pu,qu);

        % assembles new M and U matrices
        M = [pm,qm,rm];
        U = [pu,qu,ru];

        % computes DCM
        A = M*inv(U);
        
    % DCM for 3 measurements
    elseif N == 3
        A = M*inv(U);
        
    % DCM for more than 3 measurements
    else
        A = M*U'*inv(V*V');
        
    end      

end