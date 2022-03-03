%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class Definition for a GPS Constellation object
classdef GPS_Constellation < handle
    
    properties
        gps_ephem_filepath
        eph_data
        MJD_0
        GPS_MJDs
        GPS_times
        ECEFs
        ECIs
        closest4
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = GPS_Constellation(simdata)
            
            obj.gps_ephem_filepath = 'ephem.csv';
            obj.eph_data = obj.get_CSV_data(obj.gps_ephem_filepath);
            
            obj.MJD_0 = 44244; % MJD of 01/06/1980 midnight (start of GPS time)
            
            obj.GPS_MJDs = simdata.MJD_GPS;
            obj.GPS_times = obj.get_GPS_times(obj.GPS_MJDs);
            
            obj.ECEFs = 0 % obj.init_ECEFs();
            obj.ECIs = 0; % obj.init_ECIs(simdata);
            obj.closest4 = 0; % obj.init_4_closest(simdata);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get 4 pseudoranges at a specified timestep
        % 
        % Inputs: MJD - (scalar) one timestep of 'simdata.MJD_GPS'
        %         ecef_rcv - (3x1) one timestep of 'simdata.r_ecef'
        %         bias_rcv - (scalar) current estimated receiver clock bias [meters]
        %         NOISE_FLAG - (0 or 1) Add noise to measurements (1) OR
        %                               return exact measurements (0)
        %
        % Outputs: pseudoranges - 4x1 vector of pseudorange values from 4 closest GPS sats
        %          GPS_XYZB - 4x4 matrix containing ECEF position and clock
        %                     biases for each of the 4 GPS satellites. Each
        %                     row is a separate satellite, i.e.
        %                     GPS_XYZB(1,:) returns X, Y, Z, and clock bias
        %                     for the closest GPS satellite
        function [pseudoranges, GPS_XYZB] = get_pseudoranges(obj, MJD, ecef_rcv, bias_rcv, NOISE_FLAG)
            % C_LIGHT = 299792458.0;
            ind_timestep = find(MJD == obj.GPS_MJDs);
            tx_time = obj.GPS_times(ind_timestep, 2);
            SVIDs = obj.closest4(ind_timestep,:);
                        
            pseudoranges = zeros(4,1);
            GPS_XYZB = zeros(4,4);
            for ii = 1:4
                
                % Get data for current GPS sat at current timestep
                ephem = obj.eph_data(SVIDs(ii),:);
                ecef_sat = squeeze(obj.ECEFs(SVIDs(ii), ind_timestep,  1:3));
                bias_sat = obj.ECEFs(SVIDs(ii), ind_timestep, 4);
                
                % Get ionosphere stuff
                af1 = ephem.SVclockDrift;
                L1freq = 1575.42e6;
                t_rcv = tx_time - floor(tx_time/(86400*7))*86400*7;
                constant_coeffs = 1;
                Klobuchar = get_Klobuchar_coeffs(ephem, tx_time, constant_coeffs);
                I = GNSSionosphere(t_rcv,ecef_rcv,ecef_sat,Klobuchar(1,:),Klobuchar(2,:));
                
                % Calculate pseudorange 
                pseudoranges(ii) = norm(ecef_sat-ecef_rcv) + (bias_rcv-bias_sat) + I;
                
                % Also output output GPS sat's ECEF position and clock bias
                GPS_XYZB(ii,:) = [ecef_sat', bias_sat];
            end
            
            % If NOISE_FLAG set, add 0-mean Gaussian noise to measurements
            if (NOISE_FLAG == 1) 
                std_dev = 0.5;
                pseudoranges = pseudoranges + normrnd(0, std_dev, 4, 1);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initialize a 2D vector containing SVIDs of 4 closest GPS
        % satellites at each timestep
        function init_4_closest(obj, simdata)
            
            N = length(obj.GPS_times);      % Number of timesteps
            M = max(size(obj.eph_data));    % Number of GPS satellites
            closest4 = zeros(N,4);
            
            r_ecef = simdata.r_ecef;
            
            for jj = 1:N
                GPS_ECEFs = squeeze(obj.ECEFs(:, jj, 1:3));
                sc_ECEFs = repmat(r_ecef(:,jj)', 32, 1);
                
                ranges = norms(GPS_ECEFs' - sc_ECEFs')';
                [B, I] = mink(ranges, 4);
                
                closest4(jj,:) = I';
            end
            obj.closest4 = closest4;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initialize ECI positions and clock biases for all GPS satellites
        % for all timesteps. Returns a 3D array where:
        %       Dimension 1: GPS Number (SVD)
        %       Dimension 2: Timestep
        %       Dimension 3: 4-vector containing X, Y, Z, and Bias
        %
        % Note: Function takes a few seconds to run
        function init_ECIs(obj, simdata)
            
            N = length(obj.GPS_times);      % Number of timesteps
            M = max(size(obj.eph_data));    % Number of GPS satellites
            ECIs = obj.ECEFs;
            
            for jj = 1:N
                ECEF_positions = squeeze(obj.ECEFs(:,jj,1:3))';
                ROT = simdata.R_ecef2eci(:,:,jj);
                ECI_positions = ROT * ECEF_positions;
                ECIs(:,jj,1:3) = ECI_positions';
            end
            obj.ECIs = ECIs;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initialize ECEF positions and clock biases for all GPS satellites
        % for all timesteps. Returns a 3D array where:
        %       Dimension 1: GPS Number (SVD)
        %       Dimension 2: Timestep
        %       Dimension 3: 4-vector containing X, Y, Z, and Bias
        %
        % Note: Function takes a few minutes to run
        function init_ECEFs(obj)
            
            N = length(obj.GPS_times);      % Number of timesteps
            M = max(size(obj.eph_data));    % Number of GPS satellites
            ECEFs = zeros(M,N,4);
            
            for ii = 1:M
                for jj = 1:N
                    tx_time = obj.GPS_times(jj,2);
                    [ECEF_est, BIAS_est] = obj.get_sat_ECEF(obj.eph_data(ii,:), tx_time);
                    ECEFs(ii, jj, :) = [ECEF_est'; BIAS_est];
                end
            end
            obj.ECEFs = ECEFs;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Convert GPS MJD to [week, second] time format
        function GPS_times = get_GPS_times(obj, GPS_MJDs)
            GPS_MJDs_ref_1980 = GPS_MJDs - obj.MJD_0;
            GPS_total_seconds = GPS_MJDs_ref_1980 .* 86400;
            GPS_weeks = floor(GPS_total_seconds ./ (86400*7));
            GPS_secs = GPS_total_seconds - (GPS_weeks.*(86400*7));
            GPS_times = [GPS_weeks, GPS_secs];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function for calculating satellite ECEF position and clock bias 
        % given ephemeris data and transmission time
        function [ECEF_est, delta_t_SV] = get_sat_ECEF(obj, ephem, tx_time)

            % Define some constants
            MU = 3.986005e14;
            OM_E_DOT = 7.2921151467e-5;
            C = 2.99792458e8;
            F = -4.442807633e-10;
            TOL = 1e-9;

            % 19-step procedure for calculating satellite ECEF position from ephemeris
            a = ephem.sqrtA ^ 2;
            n = sqrt(MU/(a^3)) + ephem.DeltaN;
            t_k = tx_time - ephem.Toe;
            M_k = ephem.M0 + n*t_k;
            E_k = obj.MA_to_EA(M_k, ephem.Eccentricity, TOL);
            sin_v_k = sqrt(1 - ephem.Eccentricity^2)*sin(E_k) / (1 - ephem.Eccentricity*cos(E_k));
            cos_v_k = (cos(E_k)-ephem.Eccentricity) / (1 - ephem.Eccentricity*cos(E_k)); 
            v_k = atan2(sin_v_k, cos_v_k);
            phi_k = v_k + ephem.omega;
            delta_phi_k = ephem.Cus * sin(2*phi_k) + ephem.Cuc * cos(2*phi_k);
            u_k = phi_k + delta_phi_k;
            delta_r_k = ephem.Crs * sin(2*phi_k) + ephem.Crc * cos(2*phi_k);
            delta_i_k = ephem.Cis * sin(2*phi_k) + ephem.Cic * cos(2*phi_k); 
            Om_k = ephem.Omega0 - OM_E_DOT*tx_time + ephem.OmegaDot*t_k; % + -0.00000575;
            r_k = a * (1 - ephem.Eccentricity * cos(E_k)) + delta_r_k;
            i_k = ephem.Io + ephem.IDOT*t_k + delta_i_k;
            x_p = r_k * cos(u_k);
            y_p = r_k * sin(u_k); 
            X_ECEF = x_p * cos(Om_k) - y_p * cos(i_k) * sin(Om_k);
            Y_ECEF = x_p * sin(Om_k) + y_p * cos(i_k) * cos(Om_k);
            Z_ECEF = y_p * sin(i_k);
            ECEF_est = [X_ECEF, Y_ECEF, Z_ECEF];

            % Calculation of satellite's clock bias term
            delta_t_r = F * ephem.Eccentricity * ephem.sqrtA * sin(E_k);
            delta_t_SV = ephem.SVclockBias ...
                         + ephem.SVclockDrift * (tx_time - ephem.TransTime) ...
                         + ephem.SVclockDriftRate * (tx_time - ephem.TransTime)^2 ...
                         + delta_t_r ...
                         - ephem.TGD;
            delta_t_SV = delta_t_SV * C;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Helper function for reading CSV file into nice data table format
        function data = get_CSV_data(obj, filename)
            data = readtable(filename);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Helper function for iteratively solving for Eccentric Anomaly
        function EA = MA_to_EA(obj, M, e, tol)
            % Author: Joshua Geiser
            % Inputs: Mean anomaly (rad), eccentricity (--), tolerance (rad)
            % Output: Eccentric anomaly (rad)

            % Wrap M between 0 and 2 pi just to be safe
            M = wrapTo2Pi(M);

            % Safe initial guess
            E0 = M;

            % Initialize some parameters
            EA = E0;
            delta_i = 1e10;

            % Newton-Raphson method
            while (abs(delta_i) > tol)
                delta_i = -1 * (EA - e*sin(EA) - M)/(1 - e*cos(EA));
                EA = EA + delta_i;
            end
        end
        
    end 
end