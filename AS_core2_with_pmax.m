function [v_T_total,inner_products_spoofer, inner_products_jammer,y_normalized, h_m,f_m_h] = AS_core2_with_pmax(rx_signal, ssv_spoof, T, Total_snapshots, numberOfantennas,ssv_jammer, doppler_actual, actual_prns)
    % ANTI_SPOOFING Processes IQ data to detect and mitigate spoofed GPS signals
    %
    % Inputs:
    %   IQ_data           - Received IQ signal [N x numberOfantennas]
    %   ssv_spoof         - Spatial signature vector for spoofing source [4 x 1]
    %   T                 - Samples per 1ms code period (e.g., 2400)
    %   Total_snapshots   - Total number of 1ms snapshots (e.g., 30 for 30ms)
    %   numberOfantennas  - Number of antenna elements (e.g., 4)
    %
    % Outputs:
    %   v_T_total         - Processed output after null steering [N x 1]
    %   y_vectors         - Computed y vectors for each 20ms window [4 x num_windows]
    %   inner_products    - Inner products between normalized y and ssv_spoof [num_windows x 1]

    % Validate inputs
    if size(rx_signal, 2) ~= numberOfantennas
        error('Number of columns in IQ_data must match numberOfantennas.');
    end
    if length(ssv_spoof) ~= numberOfantennas
        error('Length of ssv_spoof must match numberOfantennas.');
    end

    % Initialize parameters
    N_samples = T * Total_snapshots;  % Total number of samples
    if size(rx_signal, 1) < N_samples
        error('IQ_data does not contain enough samples for specified T and Total_snapshots.');
    end
% % % % % % % % % % % % % % % % % % % % % % % % 
    
    v_T_total = zeros(N_samples, 1);  % Preallocate output


    pMaxim_snapshots = Total_snapshots;

    h_m = zeros(numberOfantennas,1);

    q_m = zeros(numberOfantennas,1);
    
    % h_m = zeros(4, 1) + 1j *zeros(4, 1);
    
    v_T_total = zeros(1, T*Total_snapshots);
    
  
% % % % % % % % % % % % % % % % % % % % % % % 
  
 %% Anti-Spoofing Processing
IQ_data = rx_signal;  % Example IQ data
 

v_T_PMax_all_prn =zeros(1,length(IQ_data));

IQ_transposed = IQ_data.';

% Initialize result storage for each antenna
theta_results_per_antenna = zeros(4, 1);  % 4 antennas and time periods

theta_results_per_antenna(1,1) = 1;


% Loop over antennas
for antenna_idx = 2:4
    theta_beta_sum_of_Kperiods_result = 0;
    Tsamples_sum_result = 0;  % Initialize sum for this particular period and antenna
    % Loop over samples in one period (T to 2*T)
    for sample_number = T:2*T
        singleIQ_result = IQ_data(sample_number, antenna_idx) * conj(IQ_data(sample_number, 1));
        Tsamples_sum_result = Tsamples_sum_result + singleIQ_result;
    end
    theta_results_per_antenna(antenna_idx, 1) = Tsamples_sum_result;
end


% Calculate phase of the complex values
phases = angle(theta_results_per_antenna(1:4, 1));


%% Beta calculation
beta_results_per_antenna = zeros(4, 1);  % 4 antennas

% Initialize Tsamples_sum_result before the loop
for antenna_idx1 = 1:4
    Tsamples_sum_result = 0;  % Initialize sum for this antenna
    for sample_number1 = T+1:2*T
        beta_singleIQ_result = IQ_data(sample_number1, antenna_idx1) * conj(IQ_data(sample_number1-T, antenna_idx1));
        Tsamples_sum_result = Tsamples_sum_result + beta_singleIQ_result;
    end
    beta_results_per_antenna(antenna_idx1, 1) = sqrt(abs(Tsamples_sum_result)); % Use abs to ensure real sqrt
end


%% Compute y = beta * e^(j*theta)
beta = beta_results_per_antenna;
theta = phases;

y = ones(4, 1);  % For i=1, e^{j theta_1} = 1
for i = 1:4
    y(i) = beta(i) * exp(1j * theta(i));
end

y_normalized = y / norm(y);


 ssv_spoof_normalized = ssv_spoof / norm(ssv_spoof);
 inner_products_spoofer= abs(y_normalized' * ssv_spoof_normalized);

    
ssv_jammer_normalized = ssv_jammer / norm(ssv_jammer);
inner_products_jammer = abs(y_normalized' * ssv_jammer_normalized);

            
% ssv_spoof_normalized = ssv_spoof / norm(ssv_spoof);

inner_product = abs(y_normalized' * ssv_spoof_normalized);

%% Null Steering Unit
y_H = y'; % Conjugate transpose of y
P_perp = eye(4) - y * (inv(y_H * y) * y_H);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % %% % % % %   POWWER MAXIM>    % % % % % % % % % % % % % %

        for Prn_num_to_Pmax = 1:length(actual_prns)

          


         for snapshot_num= 1:pMaxim_snapshots-2  
     

            P_Maxim_Tsamples_sum_result =0;
            for sample_number1 = T:2*T 

                 P_Maxim_singleIQ_result = P_perp*(IQ_data((T*snapshot_num+sample_number1), :).')* conj(IQ_data(sample_number1, antenna_idx1));

                 P_Maxim_Tsamples_sum_result = P_Maxim_Tsamples_sum_result + P_Maxim_singleIQ_result;
            end

            q_m = q_m + P_Maxim_Tsamples_sum_result*exp(-1j*2*pi*doppler_actual(Prn_num_to_Pmax)*snapshot_num);
         end

       h_m = q_m/norm(q_m);
       h_m_H = h_m';
       P_perp_H = P_perp';

       f_m_h = h_m_H*P_perp_H;
       vm =  f_m_h* IQ_transposed;
       v_T_PMax_all_prn = ( v_T_PMax_all_prn + vm);

        end

        v_T_total = v_T_PMax_all_prn.';