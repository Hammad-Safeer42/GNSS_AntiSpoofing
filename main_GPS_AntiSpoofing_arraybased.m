clc;
clear;

addpath('Library\');

% ######## Input Parameters  ############
% #######################################


Fs = 2.4e6;                  % Sampling frequency: 2.4 MHz
T_sim_per_segment = 0.440;   % Simulation duration per segment: 110 ms
T_sim = T_sim_per_segment * 3;
N_per_segment = round(T_sim_per_segment * Fs);
T = 2400;                   % One code period (1 ms at 2.4 MHz)
K = 2;                      % For 30ms
Total_snapshots = (T_sim * 1000);

numberOfantennas = 4;
fc = 1575.42e6;              % GPS L1 carrier frequency
c = 3e8;                     % Speed of light (m/s)
lambda = c / fc;             % Wavelength
d = lambda / 2;              % Element spacing (half wavelength)
acq_snapshots_factor_pmax = 5;

v_T_total = zeros(T * Total_snapshots, 1);

% 2x2 array positions
elem_pos = [0 0 0; d 0 0; 0 d 0; d d 0];


% % % % % % % % % %  PRNS and other sig Parameters   % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

actual_prns = 1:6;
spoofed_prns = [1, 2, 3, 4, 5, 6, 14, 12];

% Jammer properties
doa_az_jammer = -30;
doa_el_jammer = -40;
P_jammer_dBm = -99900;

% Power settings
P_actual_dBm = -125;
P_spoof_dBm = -122;
P_noise_dBmMHz = -115;
B_Hz = 2e6;

% Optional: Process in 11ms chunks
acq_samples = 26400;


% Define parameters for each segment
%  three segments initialized so data can be generated for three different
%  scenarios in same signal

segments = struct();
segments(1).doa_az_actual = [10, 250, 100, 160, 60, 190];
segments(1).doa_el_actual = [60, 20, 65, 30, 50, 70];
segments(1).doppler_actual = [4000, -1500, -2370, 670, 940, 110];
segments(1).doa_az_spoof = -45;
segments(1).doa_el_spoof = -65;
segments(1).doppler_spoof = [5000, 1130, 500, -800, -4000, -600, 2240, -1450];

segments(2).doa_az_actual =  [10, 250, 100, 160, 60, 190];
segments(2).doa_el_actual = [60, 20, 65, 30, 50, 70];
 segments(2).doppler_actual = [4000, -1500, -2370, 670, 940, 110];

segments(2).doa_az_spoof = -45;
segments(2).doa_el_spoof = -65;
segments(2).doppler_spoof = [5100, 1230, 600, -700, -3900, -500, 2340, -1350];

segments(3).doa_az_actual = [10, 250, 100, 160, 60, 190];
segments(3).doa_el_actual = [60, 20, 65, 30, 50, 70];
segments(3).doppler_actual = [4000, -1500, -2370, 670, 940, 110];
segments(3).doa_az_spoof = -45;
segments(3).doa_el_spoof = -65;
segments(3).doppler_spoof = [5200, 1330, 700, -600, -3800, -400, 2440, -1250];

% Initialize concatenated signals
rx_signal_total = [];
ssv_spoof_total = [];
ssv_jammer_total = [];
code_seq_total = [];


% Generate and concatenate signals


for seg = 1:3
    [rx_signal, ssv_spoof, ssv_jammer, code_seq] = simulate_gps_signal_noActualData(...
        Fs, T_sim_per_segment, fc, c, elem_pos, actual_prns, ...
        segments(seg).doa_az_actual, segments(seg).doa_el_actual, segments(seg).doppler_actual, ...
        spoofed_prns, segments(seg).doa_az_spoof, segments(seg).doa_el_spoof, segments(seg).doppler_spoof, ...
        doa_az_jammer, doa_el_jammer, P_actual_dBm, P_spoof_dBm, P_jammer_dBm, P_noise_dBmMHz, B_Hz);
    
    rx_signal_total = [rx_signal_total; rx_signal];
    ssv_spoof_total = [ssv_spoof_total; ssv_spoof];
    ssv_jammer_total = [ssv_jammer_total; ssv_jammer];
    code_seq_total = [code_seq_total; code_seq];



end




chunk_samples = 26400*acq_snapshots_factor_pmax;
T_sim_total = T_sim_per_segment * 3;
num_chunks = floor(T_sim_total * 1000 / (acq_snapshots_factor_pmax*11));



for k = 1:num_chunks
    start_idx = (k-1) * chunk_samples + 1;
    end_idx = k * chunk_samples;
    rx_signal_chunk = rx_signal_total(start_idx:end_idx, :);


% considering 10 chunks per segment, if less then need to change the
% divider

    segment_index = ceil(k /(acq_snapshots_factor_pmax*11) );   

      % Anti-spoofing processing

    Total_snapshots = acq_snapshots_factor_pmax*11;  % 11 snapshots for 11ms (11 * 2400 = 26,400 samples)

    doa_az_actual = segments(segment_index).doa_az_actual;  % Vector of 6 satellite azimuths
    doa_el_actual = segments(segment_index).doa_el_actual;  % Vector of 6 satellite elevations
    doa_az_spoof = segments(segment_index).doa_az_spoof;    % Scalar spoofer azimuth
    doa_el_spoof = segments(segment_index).doa_el_spoof;    % Scalar spoofer elevation
    % doa_az_jammer = doa_az_jammer;                          % Constant jammer azimuth
    % doa_el_jammer = doa_el_jammer;                        
   
    doa_az_all = [doa_az_actual, doa_az_spoof, doa_az_jammer];
    doa_el_all = [doa_el_actual, doa_el_spoof, doa_el_jammer];
    
    % Convert angles to radians
    doa_az_rad = deg2rad(doa_az_all);
    doa_el_rad = deg2rad(doa_el_all);
    
    x = cos(doa_el_rad) .* sin(doa_az_rad);
    y = cos(doa_el_rad) .* cos(doa_az_rad);
    z = sin(doa_el_rad);

   % if (k==1)
   dopplerss_actual_p= segments(1).doppler_actual;
   actuals_prns_p = actual_prns;
   % end

% Passing data to Anti spoofing core1 

    [v_T_total, inner_products_spoofer, inner_products_jammer, y_normalized, h, f_H] = ...
        AS_Core1_withoutPMax(rx_signal_chunk, ssv_spoof, T, Total_snapshots, numberOfantennas, ssv_jammer);

% Passing data to Anti spoofing core2 

    [v_T_total_p, inner_products_spoofer_p, inner_products_jammer_p, y_normalized_pmax, h_m, f_m_h] = ...
        AS_core2_with_pmax(rx_signal_chunk, ssv_spoof, T, Total_snapshots, numberOfantennas, ssv_jammer,dopplerss_actual_p, actuals_prns_p);
    
   

 % Acquisition on anti-spoofed signal with P maxx

   acqObj_pmax = GPSAquisition();
    acqObj_pmax.samplingFreq = Fs;
    acqObj_pmax.IF = 0;
    samples_after_p = v_T_total_p(1:acq_samples);
    [satPeaks_after_p, acqThreshold_after_p] = acqObj_pmax.step(samples_after_p);
    allDopplers_afterAS_p = acqObj_pmax.getDoppler(samples_after_p);
    detectedPRNs_after_p = find(satPeaks_after_p > acqThreshold_after_p);

     % Suppose allDopplers_afterAS is your full vector:
dopplerss_actual_p = allDopplers_afterAS_p;

% Create a logical mask for entries that are neither zero nor NaN
mask = (dopplerss_actual_p ~= 0) & ~isnan(dopplerss_actual_p);

% Apply the mask to keep only the desired values
dopplerss_actual_p = dopplerss_actual_p(mask).';



actuals_prns_p = detectedPRNs_after_p.';





% Acquisition on anti spoofed signal

    acqObj = GPSAquisition();
    acqObj.samplingFreq = Fs;
    acqObj.IF = 0;
    samples_after = v_T_total(1:acq_samples);
    [satPeaks_after, acqThreshold_after] = acqObj.step(samples_after);
    allDopplers_afterAS = acqObj.getDoppler(samples_after);
    detectedPRNs_after = find(satPeaks_after > acqThreshold_after);

     % Suppose allDopplers_afterAS is your full vector:
dopplerss_actual = allDopplers_afterAS;

% Create a logical mask for entries that are neither zero nor NaN
mask = (dopplerss_actual ~= 0) & ~isnan(dopplerss_actual);

% Apply the mask to keep only the desired values
dopplerss_actual = dopplerss_actual(mask).';



actuals_prns = detectedPRNs_after.';



% % % 

    % Acquisition on original signal (first antenna)
    acqObj_original = GPSAquisition();
    acqObj_original.samplingFreq = Fs;
    acqObj_original.IF = 0;
    samples_before = rx_signal_chunk(1:acq_samples, 1);
    [satPeaks_before, acqThreshold_before] = acqObj_original.step(samples_before);
    allDopplers_beforeAS = acqObj_original.getDoppler(samples_before);
    detectedPRNs_before = find(satPeaks_before > acqThreshold_before);






    % Beam pattern calculation
    w0 = f_m_h.';  % Beamforming weights
    steer = @(az, el) exp(-1j * 2 * pi / lambda * (elem_pos * [cosd(el)*cosd(az); cosd(el)*sind(az); sind(el)]));
    az_grid = -180:2:180;
    el_grid = 0:2:90;
    [AZ, EL] = meshgrid(az_grid, el_grid);
    r = 90 - EL;
    Gain_dB = zeros(size(AZ));
    for ii = 1:numel(AZ)
        a = steer(AZ(ii), EL(ii));
        g = w0' * a;
        Gain_dB(ii) = 20 * log10(abs(g));
    end
    Gain_dB = Gain_dB - max(Gain_dB(:));
    Gain_dB(Gain_dB < -60) = -60;
    X = r .* cosd(AZ);
    Y = r .* sind(AZ);

    % DOA estimation
    [az_sp, el_sp] = estDOAfromSSV(ssv_spoof, elem_pos, lambda);
    [az_jm, el_jm] = estDOAfromSSV(ssv_jammer, elem_pos, lambda);
    [az_y, el_y] = estDOAfromSSV(y_normalized, elem_pos, lambda);
    [az_h, el_h] = estDOAfromSSV(h_m, elem_pos, lambda);
    u_sp = [cosd(el_sp)*sind(az_sp), cosd(el_sp)*cosd(az_sp), sind(el_sp)];
    u_jm = [cosd(el_jm)*sind(az_jm), cosd(el_jm)*cosd(az_jm), sind(el_jm)];
    u_y = [cosd(el_y)*sind(az_y), cosd(el_y)*cosd(az_y), sind(el_y)];
    u_h = [cosd(el_h)*sind(az_h), cosd(el_h)*cosd(az_h), sind(el_h)];

    % Prepare data structure for visualization
    data = struct();
    data.x = x;
    data.y = y;
    data.z = z;
    data.actual_prns = actual_prns;
    data.satPeaks_before = satPeaks_before;
    data.acqThreshold_before = acqThreshold_before;
    data.satPeaks_after = satPeaks_after;
    data.acqThreshold_after = acqThreshold_after;
    data.satPeaks_after_p = satPeaks_after_p;
    data.acqThreshold_after_p = acqThreshold_after_p;

    data.X = X;
    data.Y = Y;
    data.Gain_dB = Gain_dB;
    data.doa_el_actual = abs(doa_el_actual);
    data.doa_az_actual = doa_az_actual;
    data.doa_el_spoof = abs(doa_el_spoof);
    data.doa_az_spoof = doa_az_spoof;
    data.u_sp = u_sp;
    data.u_jm = u_jm;
    data.u_y = u_y;
    data.u_h = u_h;
    data.chunk = k; 



    % Update visualization
    plot_all_data(data);
    pause(0.5);  

end

function [best_az, best_el] = estDOAfromSSV(ssv, elem_pos, lambda)
    az_grid = -180:1:180;
    el_grid = -90:1:90;
    best_val = -inf;
    best_az = 0; best_el = 0;
    for az = az_grid
        for el = el_grid
            u = [cosd(el)*cosd(az); cosd(el)*sind(az); sind(el)];
            a = exp(1j*2*pi*(elem_pos * u)/lambda);
            m = abs(a' * ssv);
            if m > best_val
                best_val = m;
                best_az = az;
                best_el = el;
            end
        end
    end
end


