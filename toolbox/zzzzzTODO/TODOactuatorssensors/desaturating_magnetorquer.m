%==========================================================================
%
% desaturating_magnetorquer  Simulates the effect of a desaturating
% magnetorquer.
%
%   w = desaturating_magnetorquer(w,w_max)
%
% Copyright (c) 2021 Luke Neise
% Last Update: 2021-06-02
%
%--------------------------------------------------------------------------
%
% -------
% INPUTS:
% -------
%   w       (nwx1) [rad/s] angular velocity of each gyroscopic actuator
%   w_max   (nwx1) [rad/s] max angular velocity of each gyroscopic actuator
%
% --------
% OUTPUTS:
% --------
%   w       (nwx1) [rad/s] angular velocity of each gyroscopic actuator
%                          after desaturation has been accounted for
%
%==========================================================================
function w = desaturating_magnetorquer(w,w_max)

    % convert current rotational velocity of reaction wheels to rpm
    w_rpm = w*(60/(2*pi));

    % check against saturation limit and desaturate if exceeded
    for j = 1:length(w)
        if w_rpm(j) >= w_max
            w(j) = 0;
        end
    end

end