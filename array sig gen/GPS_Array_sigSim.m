clc;
clear;



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




% % % % % % % % % % % % % % 
    %  " rx_signal "   is your desired simulated signal of gps (no nav data
    %  ) on 2x2 array with different sat signals , spoofer and jammer
% % % % % % % % % % % % % 
