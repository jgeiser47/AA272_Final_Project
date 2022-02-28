%==========================================================================
%
% true_state  Ground truth state.
%
%   x_true = true_state(simdata)
%
% Author: Tamas Kis
% Last Update: 2022-02-21
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   simdata - (1×1 struct) simulation data
%
% -------
% OUTPUT:
% -------
%   x_true  - (6×N double) ground truth state
%               --> 1-3. position resolved in ECI frame [m]
%               --> 4-6. inertial velocity resolved in ECI frame [m/s]
%
% -----
% NOTE:
% -----
%   --> N = number of iterations (i.e. length of time vector)
%
%==========================================================================
function x_true = true_state(simdata)
    x_true = [simdata.r_eci;
              simdata.v_eci];
end