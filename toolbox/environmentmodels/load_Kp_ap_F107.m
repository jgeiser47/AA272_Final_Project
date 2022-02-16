%==========================================================================
%
% load_Kp_ap_F107  Loads space weather data.
%
%   Kp_ap_F107 = load_Kp_ap_F107
%   Kp_ap_F107 = load_Kp_ap_F107(MJD_UTC0,duration)
%
% Author: Tamas Kis
% Last Update: 2022-02-15
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   MJD_UTC0    - (1×1 double) (OPTIONAL) initial UTC [MJD]
%   duration    - (1×1 double) (OPTIONAL) simulation duration [h]
%
% -------
% OUTPUT:
% -------
%   Kp_ap_F107  - (N×20 double) space weather data; columns:
%                    1. MJD_UT1  - UT1 (Universal Time 1)[MJD]
%                    2. Kp1      - planetary index (00:00-03:00) [-]
%                    3. Kp2      - planetary index (03:00-06:00) [-]
%                    4. Kp3      - planetary index (06:00-09:00) [-]
%                    5. Kp4      - planetary index (09:00-12:00) [-]
%                    6. Kp5      - planetary index (12:00-15:00) [-]
%                    7. Kp6      - planetary index (15:00-18:00) [-]
%                    8. Kp7      - planetary index (18:00-21:00) [-]
%                    9. Kp8      - planetary index (21:00-00:00) [-]
%                   10. ap1      - planetary amplitude (00:00-03:00) [γ]
%                   11. ap2      - planetary amplitude (03:00-06:00) [γ]
%                   12. ap3      - planetary amplitude (06:00-09:00) [γ]
%                   13. ap4      - planetary amplitude (09:00-12:00) [γ]
%                   14. ap5      - planetary amplitude (12:00-15:00) [γ]
%                   15. ap6      - planetary amplitude (15:00-18:00) [γ]
%                   16. ap7      - planetary amplitude (18:00-21:00) [γ]
%                   17. ap8      - planetary amplitude (21:00-00:00) [γ]
%                   18. Ap       - daily planetary amplitude [γ]
%                   19. F107     - 10.7 cm solar flux for prev. day [SFU]
%                   20. F107_avg - centered 81-day average of F10.7 [SFU]
%
% -----
% NOTE:
% -----
%   --> N = # of data entries
%
%==========================================================================
function Kp_ap_F107 = load_Kp_ap_F107(MJD_UTC0,duration)
    
    % loads full data set
    Kp_ap_F107 = struct2array(load('Kp_ap_F107.mat'));
    
    % trims data set to only keep relevant data for simulation
    if nargin ~= 0

        % start and end dates of data to keep (from 3 days before start to 
        % end of simulation)
        MJD_start = floor(MJD_UTC0)-3;
        MJD_end = ceil(MJD_UTC0+duration/24);
    
        % start and end indices
        i_start = interval_search(Kp_ap_F107(:,1),MJD_start,false);
        [~,i_end] = interval_search(Kp_ap_F107(:,1),MJD_end,false);
    
        % removes unnecessary data
        Kp_ap_F107 = Kp_ap_F107(i_start:i_end,:);

    end

end