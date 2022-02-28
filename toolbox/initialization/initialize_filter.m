%==========================================================================
%
% initialize_filter  Filter initialization.
%
%   filter = initialize_filter(UTC_start,duration,sim_start,N)
%
% Author: Tamas Kis
% Last Update: 2022-02-21
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   UTC_start   - (1×6 double) initial UTC [y,mo,d,h,m,s]
%   duration    - (1×1 double) simulation duration [hr]
%   sim_start   - (OPTIONAL) (1×1 double) initial simulation time [s]
%   N           - (OPTIONAL) (1×1 double) maximum degree/order of gravity
%                  model to use
%
% -------
% OUTPUT:
% -------
%   filter      - (1×1 struct) structure to store all filter settings
%       • data - (1×1 struct) data for orbit propagation/determination
%           ‣ C         - (1×1 function_handle) gravitational coefficients
%           ‣ S         - (1×1 function_handle) gravitational coefficients
%           ‣ eop       - (1×1 struct) Earth orientation parameters
%           ‣ XYs_iau06 - (1×1 struct) additional data for calculating
%                         Earth orientation
%           ‣ sw1       - (1×1 struct) 1st half of space weather data
%           ‣ sw2       - (1×1 struct) 2nd half of space weather data
%           ‣ nrlm_data - (1×1 struct) data needed for NRLMSISE-00
%                         atmospheric model
%     	• t0 - (1×1 struct) initial times
%           ‣ t0       - (1×1 double) simulation time [s]
%           ‣ MJD_GPS0 - (1×1 double) GPS time [MJD]
%           ‣ MJD_TAI0 - (1×1 double) TAI (International Atomic Time) [MJD]
%           ‣ MJD_TT0  - (1×1 double) TT (Terrestrial Time) [MJD]
%           ‣ MJD_UT10 - (1×1 double) UT1 (Universal Time 1) [MJD]
%           ‣ MJD_UTC0 - (1×1 double) UTC (Universal Coordinated Time) 
%                        [MJD]
%
%==========================================================================
function filter = initialize_filter(UTC_start,duration,sim_start,N)

    % -------------------------------------------------------
    % Defaults optional inputs to empty vectors if not input.
    % -------------------------------------------------------

    if (nargin < 2), duration = []; end
    if (nargin < 3), sim_start = []; end
    if (nargin < 4), N = []; end
    
    % --------------
    % Initial times.
    % --------------

    % defaults initial simulation time to 0 if not input
    if (nargin < 6) || isempty(sim_start)
        sim_start = 0;
    end
    
    % initialize structure of initial times
    t0 = initialize_time(UTC_start,sim_start);

    % -----
    % Data.
    % -----

    data = initialize_data(t0.MJD_UTC0,duration,N);
    
    % -------------------------------------
    % Package overall propagator structure.
    % -------------------------------------

    filter.data = data;
    filter.t0 = t0;

end