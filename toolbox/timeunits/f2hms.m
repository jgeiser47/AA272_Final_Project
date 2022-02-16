%==========================================================================
%
% f2hms  Hours, minutes, and seconds from fraction of day.
%
%   [hh,mm,ss] = f2hms(f_DD)
%
% See also cal2doy, cal2mjd, hms2f, jd2mjd, jd2t, mjd2cal, mjd2jd, mjd2t.
%
% Author: Tamas Kis
% Last Update: 2022-02-01
%
% REFERENCES:
%   [1] Vallado, "Fundamentals of Astrodynamics and Applications", 4th Ed.,
%       (p. 199)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f_DD    - (1×1 double) fraction of day [d]
%
% -------
% OUTPUT:
% -------
%   hh      - (1×1 double) hours [h]
%   mm      - (1×1 double) minutes [m]
%   ss      - (1×1 double) seconds [s]
%
%==========================================================================
function [hh,mm,ss] = f2hms(f_DD)
    
    % hours
    hh = floor(24*f_DD);

    % minutes
    mm = floor(60*(24*f_DD-hh));

    % seconds
    ss = 3600*(24*f_DD-hh-(mm/60));
    
end