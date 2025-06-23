function [rx_signal, ssv_spoof, ssv_jammer, code_seq] = simulate_gps_signal_noActualData(Fs, T_sim, fc, c, elem_pos, actual_prns, doa_az_actual, doa_el_actual, doppler_actual, spoofed_prns, doa_az_spoof, doa_el_spoof, doppler_spoof, doa_az_jammer, doa_el_jammer, P_actual_dBm, P_spoof_dBm, P_jammer_dBm, P_noise_dBmMHz, B_Hz)


% SIMULATE_GPS_SIGNAL_NOACTUALDATA Simulates GPS signals with actual, spoofed PRNs, and a jammer
%
% Inputs:
%   Fs                    - Sampling frequency (Hz)
%   T_sim                 - Simulation duration (seconds)
%   fc                    - Carrier frequency (Hz)
%   c                     - Speed of light (m/s)
%   elem_pos              - Antenna element positions (m) [numberOfantennas x 3]
%   actual_prns           - Actual PRN numbers (vector)
%   doa_az_actual         - Azimuth angles for actual PRNs (degrees)
%   doa_el_actual         - Elevation angles for actual PRNs (degrees)
%   doppler_actual        - Doppler shifts for actual PRNs (Hz)
%   spoofed_prns          - Spoofed PRN numbers (vector)
%   doa_az_spoof          - Common azimuth angle for spoofed PRNs (degrees)
%   doa_el_spoof          - Common elevation angle for spoofed PRNs (degrees)
%   doppler_spoof         - Doppler shifts for spoofed PRNs (Hz)
%   doa_az_jammer         - Jammer azimuth angle (degrees)
%   doa_el_jammer         - Jammer elevation angle (degrees)
%   P_actual_dBm          - Power per actual signal in dBm
%   P_spoof_dBm           - Power per spoofed signal in dBm
%   P_jammer_dBm          - Jammer power in dBm
%   P_noise_dBmMHz        - Noise power spectral density in dBm/MHz
%   B_Hz                  - Bandwidth in Hz
%
% Outputs:
%   rx_signal             - Simulated received signal [N x numberOfantennas]
%   ssv_spoof             - Spatial signature vector for spoofing source [4 x 1]
%   ssv_jammer            - Spatial signature vector for jammer [4 x 1]
%   code_seq              - C/A code sequences [15 x 1023]



% Derived parameters
N = round(T_sim * Fs);              % Number of samples
lambda = c / fc;                    % Wavelength
numberOfantennas = size(elem_pos, 1); % Number of antenna elements

% Power control: Convert dBm to watts
B_MHz = B_Hz / 1e6;                 % Bandwidth in MHz
P_noise_dBm = P_noise_dBmMHz + 10 * log10(B_MHz); % Total noise power in dBm
P_actual = 10^((P_actual_dBm - 30)/10); % dBm to watts
P_spoof = 10^((P_spoof_dBm - 30)/10);
P_jammer = 10^((P_jammer_dBm - 30)/10);
P_noise = 10^((P_noise_dBm - 30)/10);


% Set signal amplitudes (amplitude = sqrt(power))
signal_ampl = [sqrt(P_actual) * ones(1, numel(actual_prns)), sqrt(P_spoof) * ones(1, numel(spoofed_prns))];

% Generate C/A codes for PRNs 1 to 15
num_chips = 1023;                   % C/A code length
num_corr_prns = 15;                 % PRNs to generate
code_seq = zeros(num_corr_prns, num_chips);
g2_taps = {
    [2 6], [3 7], [4 8], [5 9], [1 9], [2 10], [1 8], ...
    [2 9], [3 10], [2 3], [3 4], [5 6], [6 7], [7 8], ...
    [8 9]
};

for prn = 1:num_corr_prns
    G1 = ones(1, 10);
    G2 = ones(1, 10);
    taps = g2_taps{prn};
    code_bits = zeros(1, num_chips);
    for i = 1:num_chips
        g1_out = G1(10);
        g2_out = xor(G2(taps(1)), G2(taps(2)));
        code_bits(i) = xor(g1_out, g2_out);
        new_bit_G1 = xor(G1(3), G1(10));
        new_bit_G2 = xor(xor(xor(G2(2), G2(3)), G2(6)), xor(xor(G2(8), G2(9)), G2(10)));
        G1 = [new_bit_G1, G1(1:9)];
        G2 = [new_bit_G2, G2(1:9)];
    end
    code_seq(prn, :) = 2 * code_bits - 1;  % Map 0/1 to ±1
end

% Prepare signal properties
signal_prn = [actual_prns, spoofed_prns];
signal_az = [doa_az_actual, doa_az_spoof * ones(1, length(spoofed_prns))];
signal_el = [doa_el_actual, doa_el_spoof * ones(1, length(spoofed_prns))];
signal_doppler = [doppler_actual, doppler_spoof];
signal_freqs = 0 - signal_doppler; % f_IF = 0
% signal_freqs = 0 + signal_doppler; % f_IF = 0
% Compute SSV for spoofing source
az_rad_spoof = deg2rad(doa_az_spoof);
el_rad_spoof = deg2rad(doa_el_spoof);
dir_vec_spoof = [cos(el_rad_spoof)*cos(az_rad_spoof), ...
                 cos(el_rad_spoof)*sin(az_rad_spoof), ...
                 sin(el_rad_spoof)];
ssv_spoof = zeros(4, 1);
for elem = 1:4
    proj = dot(elem_pos(elem, :), dir_vec_spoof);
    phase_offset = 2 * pi * proj / lambda;
    ssv_spoof(elem) = exp(1j * phase_offset);
end

% Compute SSV for jammer
az_rad_jammer = deg2rad(doa_az_jammer);
el_rad_jammer = deg2rad(doa_el_jammer);
dir_vec_jammer = [cos(el_rad_jammer)*cos(az_rad_jammer), ...
                  cos(el_rad_jammer)*sin(az_rad_jammer), ...
                  sin(el_rad_jammer)];
ssv_jammer = zeros(4, 1);
for elem = 1:4
    proj = dot(elem_pos(elem, :), dir_vec_jammer);
    phase_offset = 2 * pi * proj / lambda;
    ssv_jammer(elem) = exp(1j * phase_offset);
end

% Simulate received signal
t = (0:N-1).' / Fs;                % Time vector
rx_signal = zeros(N, numberOfantennas);  % Initialize as complex
code_rate = 1.023e6;               % C/A code chip rate
chip_idx = mod(floor(t * code_rate), num_chips) + 1;

% Add actual and spoofed signals
for sig_idx = 1:length(signal_prn)
    prn = signal_prn(sig_idx);
    code_wave = code_seq(prn, chip_idx).';  % N x 1 code waveform
    f_sig = signal_freqs(sig_idx);
    az_rad = deg2rad(signal_az(sig_idx));
    el_rad = deg2rad(signal_el(sig_idx));
    dir_vec = [cos(el_rad)*cos(az_rad), cos(el_rad)*sin(az_rad), sin(el_rad)];
    for elem = 1:numberOfantennas
        proj = dot(elem_pos(elem, :), dir_vec);
        phase_offset = 2 * pi * proj / lambda;
        phi_t = 2 * pi * f_sig * t + phase_offset;
        sig_wave_elem = signal_ampl(sig_idx) * code_wave .* exp(1j * phi_t);
        rx_signal(:, elem) = rx_signal(:, elem) + sig_wave_elem;
    end
end

% Add jammer signal (broadband Gaussian noise)
jammer_ampl = sqrt(P_jammer / 2); % Split power between I and Q
for elem = 1:numberOfantennas
    proj = dot(elem_pos(elem, :), dir_vec_jammer);
    phase_offset = 2 * pi * proj / lambda;
    jammer_noise = jammer_ampl * (randn(N, 1) + 1j * randn(N, 1)) * exp(1j * phase_offset);
    rx_signal(:, elem) = rx_signal(:, elem) + jammer_noise;
end

% Add thermal noise
for elem = 1:numberOfantennas
    noise = sqrt(P_noise / 2) * (randn(N, 1) + 1j * randn(N, 1));
    rx_signal(:, elem) = rx_signal(:, elem) + noise;
end
end