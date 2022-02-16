%==========================================================================
%
% space_weather  Space weather data for the Jacchia-Roberts and NRLMSISE-00
% atmospheric models.
%
%   sw = space_weather(MJD_UT1,Kp_ap_F017)
%
% Author: Tamas Kis
% Last Update: 2022-02-15
%
% REFERENCES:
%   [1] Vallado, "Fundamentals of Astrodynamics and Applications", 4th Ed.
%       (pp. 560, 604, 1001-1002, 1053-1055)
%   [2] Long et al., "Goddard Trajectory Determination System (GTDS)
%       Mathematical Theory: Revision 1" (pp. 4-37 to 4-38)
%   [3] https://celestrak.com/SpaceData/
%   [4] https://celestrak.com/SpaceData/SpaceWx-format.php
%   [5] nrlmsise-00.h 
%       (https://git.linta.de/?p=~brodo/nrlmsise-00.git;a=summary)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   MJD_UT1     - (1×1 double) UT1 (Universal Time 1) [MJD]
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
%                   19. F107     - daily 10.7 cm solar flux [SFU]
%                   20. F107_avg - centered 81-day average of F10.7 [SFU]
%
% -------
% OUTPUT:
% -------
%   sw          - (1×1 struct) space weather data for Jacchia-Roberts and 
%                 NRLMSISE-00
%       • Ap        - (1×1 double) daily planetary amplitude [γ]
%       • ap_array  - (1×7 double) planetary amplitude array for
%                     NRLMSISE-00 atmospheric model
%           ‣ 1. Ap       - daily planetary amplitude [γ]
%           ‣ 2. ap       - planetary amplitude (current time) [γ]
%           ‣ 3. ap3      - planetary amplitude (3 hours before) [γ]
%           ‣ 4. ap6      - planetary amplitude (6 hours before) [γ]
%           ‣ 5. ap9      - planetary amplitude (9 hours before) [γ]
%           ‣ 6. ap_12_33 - planetary amplitude (average of 12-33 hours 
%                           before) [γ]
%           ‣ 7. ap_36_57 - planetary amplitude (average of 36-57 hours 
%                           before) [γ]
%       • F107      - (1×1 double) daily 10.7 cm solar flux for previous
%                     day [SFU]
%       • F107_avg  - (1×1 double) centered 81-day average of F10.7 for
%                     previous day [SFU]
%       • Kp        - (1×1 double) planetary index 3 hours before current
%                     time [-]
%
% -----
% NOTE:
% -----
%   --> The F10.7 parameters are NOT adjusted to 1 AU.
%
%==========================================================================
function sw = space_weather(MJD_UT1,Kp_ap_F107)
    
    % extracts data from Kp_ap_F107
    MJD_UT1_data = Kp_ap_F107(:,1);
    Kp_data = Kp_ap_F107(:,2:9);
    ap_data = Kp_ap_F107(:,10:17);
    Ap_data = Kp_ap_F107(:,18);
    F107_data = Kp_ap_F107(:,19);
    F107_avg_data = Kp_ap_F107(:,20);

    % --------
    % Indices.
    % --------

    % index for current day
    today = interval_search(MJD_UT1_data,MJD_UT1,false);

    % index for previous day
    yesterday = today-1;

    % ---------------------------
    % Parameters for both models.
    % ---------------------------

    % 10.7 cm solar flux for previous day [SFU]
    F107 = F107_data(yesterday);

    % centered 81-day average of F10.7 for previous day [SFU]
    F107_avg = F107_avg_data(yesterday);

    % -------------------------------------
    % Parameters for Jacchia-Roberts model.
    % -------------------------------------

    % planetary index 3 hours before current time
    Kp = t2Kp(MJD_UT1-3/24,MJD_UT1_data,Kp_data);
    
    % ---------------------------------
    % Parameters for NRLMSISE-00 model.
    % ---------------------------------

    % daily planetary amplitude [γ]
    Ap = Ap_data(today);

    % planetary amplitudes until 57 hours before current time [γ]
    ap = zeros(20,1);
    for i = 1:20
        ap(i) = t2ap(MJD_UT1-(3*(i-1))/24,MJD_UT1_data,ap_data);
    end
    
    % relevant planetary amplitudes for NRLMSISE-00 model [γ]
    ap0 = ap(1);                % current time
    ap3 = ap(2);                % 3 hours before current time
    ap6 = ap(3);                % 6 hours before current time
    ap9 = ap(4);                % 9 hours before current time
    ap_12_33 = mean(ap(5:12));  % avg. of 12-33 hours before current time
    ap_36_57 = mean(ap(13:20)); % avg. of 36-57 hours before current time

    % planetary amplitude array for NRLMSISE-00 atmospheric model [γ]
    ap_array = [Ap,ap0,ap3,ap6,ap9,ap_12_33,ap_36_57];

    % -------------------
    % Packages structure.
    % -------------------
    
    sw.Ap = Ap;
    sw.ap_array = ap_array;
    sw.F107 = F107;
    sw.F107_avg = F107_avg;
    sw.Kp = Kp;
    
end