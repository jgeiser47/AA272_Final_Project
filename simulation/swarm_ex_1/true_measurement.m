%==========================================================================
%
% true_measurement  Ground truth measurements.
%
%   y = true_measurement(x_true,simdata,prop)
%
% Author: Tamas Kis
% Last Update: 2022-02-21
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x_true  - (6×N double) ground truth state TODO
%   R       - (6×6 double) measurement noise covariance
%   simdata - (1×1 struct) simulation data TODO
%   filter  - (1×1 struct) TODO
%
% -------
% OUTPUT:
% -------
%   y       - (6×N double) ground truth measurement TODO
%
%==========================================================================
function y = true_measurement(x_true,R,simdata,filter)

    % length of time vector
    N = length(simdata.t);

    % zero vector for noise generation
    vzero = zeros(6,1);

    % sample measurement noise for each sample time
    v = zeros(6,N);
    for i = 1:N
        v(:,i) = mvnrnd(vzero,R);
    end

    % time history of ECEF state
    X_ecef = [simdata.r_ecef;
              simdata.v_ecef];

    % measurement
    y = X_ecef+v;
    
end