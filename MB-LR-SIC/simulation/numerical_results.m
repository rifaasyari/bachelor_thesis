function numerical_results(K,N_txi,scenario,var_s,modulation,runs)
% BER performance of MB-LR-SIC algorithm
% INPUT K (int): Number of users
%       N_txi (int): Number of transmit antennas for each user
%       scenario (int): 
%           1: I.i.d Gaussian random fading channels with zero mean and unit
%               variance
%           2: Realistic propagation conditions with path loss and 
%               correlated antennas
%       var_s (int): Signal power
%       modulation (string): 'QPSK' or '16-QAM'
%       runs (int): Number of simulation runs per SNR

% +++++++++++++++++++++++++ Begin configuration ++++++++++++++++++++++++++

addpath('../algorithm');

% K = 3;  % number of active user

% N_txi = 2;  % number of antennas for user i
N_t = K*N_txi;  % total number of transmit antennas
N_r = N_t;  % number of total receive antennas 

% scenario = 2;  % Choose channel scenario
% 1: Independt and identically distributed (i.i.d.) random fading channels
% 2: Realistic propagation conditions and correlated antennas

% Parameters for channel scenario 2
L_p = 0.7;  % Base power path loss
sigma = 6;  % shadowing spread in dB
correlation_tx = 0.2;  % correlation index of neighboring transmit antennas
correlation_rx = correlation_tx;  % correlation index of neighboring receive
                                  % antennas

channel_estimation = false;  % If false, assume perfectly known channel

% modulation = 'QPSK';  % 'QPSK', '16-QAM'

% var_s = 1;  % Signal power

% runs = 10;
train_symbols = 50*N_t;  % Number of symbols used for channel estimation
transmitted_symbols = 50*N_t;

snr_step = 5;

ml_detection = true;
mblrsic_dec = true;
lr_sic_dec = true;
mbsic_dec = true;
sic_dec = true;
lr_mmse_dec = false;

% +++++++++++++++++++++++++++ End configuration ++++++++++++++++++++++++++

assert(N_r >= N_t);  % Make sure to use no less receive than transmit 
                     % antennas 
                     
 start_time = datestr(datetime('now'), 'dd.mm.yyyy HH:MM:SS');
 fprintf('\nStarted at %s\n\n', start_time);
                     
if strcmp(modulation, '16-QAM')
    % var_s = 10;
    decode = @decode_16QAM;
    SNR_max = 20;
    constellation = create_16_QAM(var_s);
elseif strcmp(modulation, 'QPSK')
    % var_s = 2;
    decode = @decode_qpsk;
    SNR_max = 15;
    constellation = create_MQAM(4,var_s);  % QPSK = 4-QAM constellation
else
    error('Invalid modulation %s', modulation);
end
                     
if channel_estimation == true
    transmitted_symbols = transmitted_symbols + train_symbols;
end

R_tx = correlation_tx.^(((0:-1:-(N_t-1))'+(0:N_t-1)).^2);  % Correlation matrix
% for transmit antennas
R_rx = correlation_rx.^(((0:-1:-(N_r-1))'+(0:N_r-1)).^2);  % Correlation matrix
% for receive antennas
alpha_k = sqrt(L_p);
beta_k = 10.^(sigma*randn(N_r,N_t)/10);

if channel_estimation == true
    train_signal = constellation(mod(0:N_t-1,numel(constellation))+1)';  
    % Signal used for channel estimation 
end
                                   
BER_mblrsic = zeros(1, SNR_max/snr_step+1);
BER_mbsic = zeros(1, SNR_max/snr_step+1);  % MB-SIC detector
BER_lrsic = zeros(1, SNR_max/snr_step+1);  % LR-SIC detector
BER_sic = zeros(1, SNR_max/snr_step+1);  % SIC detector
BER_lrmmse = zeros(1, SNR_max/snr_step+1);  % LR-MMSE detector
BER_ml = zeros(1, SNR_max/snr_step+1);

for SNR = 0:snr_step:SNR_max
    
    N = N_t * var_s * 10^(-SNR/10);  % Noise power. SNR definded as 
                                     % 10*log10(N_t*var_s/N)
      
    % var_s = N/N_t * 10^(SNR/10);
                                     
    BERs_mblrsic = zeros(1, runs);
    BERs_mbsic = zeros(1, runs);
    BERs_lrsic = zeros(1,runs);
    BERs_sic = zeros(1,runs);
    BERs_lrmmse = zeros(1,runs);
    BERs_ml = zeros(1, runs);
    
    for i = 1:runs  % simulation runs
        H = 1/sqrt(2) * (randn(N_r, N_t) + 1j * randn(N_r, N_t)); 
        % channels for scenario 1
        
        if scenario == 2
            G = sqrtm(R_rx)*H*sqrtm(R_tx);
            H = alpha_k*beta_k.*G;
        end
        
        BE_mblrsic = 0;  % Bit Errors for MB-LR-SIC detector
        BE_mbsic = 0;
        BE_lrsic = 0;
        BE_sic = 0;
        BE_lrmmse = 0;
        BE_ml = 0;  % Bit Errors for ML detector
        
        H_orig = H;

        for k = 1:transmitted_symbols
            
            n = sqrt(N/2) * (randn(N_r,1)+1j*randn(N_r,1));  % noise vector
            
            if k <= train_symbols && channel_estimation == true
                y = H_orig*train_signal + n;
                
                if k == 1
                    [H, D, P] = ls_channel_estimation_mimo(y, ... 
                        train_signal, 0.75, 1e-6);
                else
                    [H, D, P] = ls_channel_estimation_mimo(y, ...
                        train_signal, 0.75, 1e-6, D, P);
                end
                
            else
            
                s = constellation(randi([1,numel(constellation)], 1, N_t))';  
                % Generate signal vector

                b = decode(s, constellation);  % Transmitted bits

                y = H * s + n;

                if mblrsic_dec == true
                    s_mblrsic = mb_lr_sic(y,H,var_s,N,modulation,constellation);
                    b_mblrsic = decode(s_mblrsic, constellation);
                    BE_mblrsic = BE_mblrsic + sum(b ~= b_mblrsic);
                end

                if lr_sic_dec == true
                    s_lrsic = lr_sic(y,H,var_s,N,modulation,constellation,true);
                    b_lrsic = decode(s_lrsic, constellation);
                    BE_lrsic = BE_lrsic + sum(b ~= b_lrsic);
                end

                if mbsic_dec == true
                    s_mbsic = mb_lr_sic(y,H,var_s,N,modulation,constellation,false);
                    b_mbsic = decode(s_mbsic, constellation);
                    BE_mbsic = BE_mbsic + sum(b ~= b_mbsic);
                end

                if sic_dec == true
                    s_sic = lr_sic(y,H,var_s,N,modulation,constellation,false);
                    b_sic = decode(s_sic, constellation);
                    BE_sic = BE_sic + sum(b ~= b_sic);
                end

                if lr_mmse_dec == true
                    s_lrmmse = lr_mmse(y,H,var_s,N,modulation,constellation,true);
                    b_lrmmse = decode(s_lrmmse, constellation);
                    BE_lrmmse = BE_lrmmse + sum(b ~= b_lrmmse);
                end

                if ml_detection == true
                    s_ml = ml_mimo(y,H,constellation);
                    b_ml = decode(s_ml, constellation);
                    BE_ml = BE_ml + sum(b ~= b_ml);
                end
            end
        end
        
        BERs_mblrsic(i) = BE_mblrsic / (transmitted_symbols * N_t * ... 
            log2(numel(constellation)));
        BERs_mbsic(i) = BE_mbsic / (transmitted_symbols * N_t * ... 
            log2(numel(constellation)));
        BERs_lrsic(i) = BE_lrsic / (transmitted_symbols * N_t * ... 
            log2(numel(constellation)));
        BERs_sic(i) = BE_sic / (transmitted_symbols * N_t * ... 
            log2(numel(constellation)));
        BERs_ml(i) = BE_ml / (transmitted_symbols * N_t * ... 
            log2(numel(constellation)));
        BERs_lrmmse(i) = BE_lrmmse / (transmitted_symbols * N_t * ... 
            log2(numel(constellation)));
        
%         if mblrsic_dec == true
%             fprintf('%d MB-LR-SIC, BER: %.4f\n', ... 
%                 i, BERs_mblrsic(i));
%         end
%         if lr_sic_dec == true
%             fprintf('%d LR-SIC, BER: %.4f\n', i, BERs_lrsic(i));
%         end
%         if mbsic_dec == true
%             fprintf('%d MB-SIC, BER: %.4f\n', i, BERs_mbsic(i));
%         end
%         if sic_dec == true
%             fprintf('%d SIC, BER: %.4f\n', i, BERs_sic(i));
%         end
%         if lr_mmse_dec == true
%             fprintf('%d LR-MMSE, BER: %.4f\n', i, BERs_lrmmse(i));
%         end
%         if ml_detection == true
%             fprintf('%d ML, BER: %.4f\n\n', i, BERs_ml(i));
%         end
        
    end
    
    fprintf('%s\n', datestr(datetime('now'), 'dd.mm.yyyy HH:MM:SS'));
    BER_mblrsic(SNR/snr_step+1) = mean(BERs_mblrsic);
    if mblrsic_dec == true
        fprintf('SNR = %d, BER MB-LR-SIC: %.4f\n', ... 
            SNR, BER_mblrsic(SNR/snr_step+1));
    end
    BER_mbsic(SNR/snr_step+1) = mean(BERs_mbsic);
    BER_lrsic(SNR/snr_step+1) = mean(BERs_lrsic);
    BER_sic(SNR/snr_step+1) = mean(BERs_sic);
    BER_lrmmse(SNR/snr_step+1) = mean(BERs_lrmmse);
    BER_ml(SNR/snr_step+1) = mean(BERs_ml);
    if ml_detection == true
        fprintf('SNR = %d, BER ML: %.4f\n\n', ... 
            SNR, BER_ml(SNR/snr_step+1));
    end
end

figure
if mblrsic_dec == true
    semilogy(0:snr_step:SNR_max, BER_mblrsic, 'r-s', ... 
        'LineWidth', 2, ...
        'DisplayName', 'MB-LR-SIC');
    hold on
end
if mbsic_dec == true
    semilogy(0:snr_step:SNR_max, BER_mbsic, 'g-o', ... 
        'LineWidth', 2, ...
        'DisplayName', 'MB-SIC');
    hold on
end
if lr_sic_dec == true
    semilogy(0:snr_step:SNR_max, BER_lrsic, 'b-p', ... 
        'LineWidth', 2, ...
        'DisplayName', 'LR-SIC');
    hold on
end
if sic_dec == true
    semilogy(0:snr_step:SNR_max, BER_sic, 'm-*', ... 
        'LineWidth', 2, ...
        'DisplayName', 'SIC');
    hold on
end
if lr_mmse_dec == true
    semilogy(0:snr_step:SNR_max, BER_lrmmse, 'b-', ... 
        'LineWidth', 2, ...
        'DisplayName', 'LR-MMSE');
    hold on
end
if ml_detection == true
    semilogy(0:snr_step:SNR_max, BER_ml, 'k-', ... 
        'LineWidth', 2, ...
        'DisplayName', 'ML');
end

title(sprintf('MU-MIMO System, $N_{txi}$=%d, $N_r$=%d, L=%d Branches and K=%d users', ...
    N_txi, N_r, N_t, K));

xlabel('SNR [dB]');
xticks(0:snr_step:SNR_max);
ylabel('BER');
ylim([0.5e-5, 1e0]);

grid on;

legend('show', 'Location', 'southwest');

annot = sprintf('Transmissions: %d\nSimulation runs: %d\nSignal power: %d\nScenario: %d', ... 
    transmitted_symbols, runs,var_s,scenario);
textbox = uicontrol('style', 'text');
set(textbox, 'String', annot);
set(textbox, 'Units', 'characters');
set(textbox, 'Position', [30, 4, 20, 4]);

savefig(sprintf('BER_%s_%s.fig', ... 
    modulation, datestr(datetime('now'), 'dd-mm-yyyy_HH-MM-SS')));
close

% filename = sprintf('BER_%s_%s.tikz', ... 
%     modulation, datestr(datetime('now'), 'dd-mm-yyyy_HH-MM-SS'));
% matlab2tikz(filename, ...
%             'height', '\fheight', 'width', '\fwidth', 'parseStrings', false);


end_time = datestr(datetime('now'), 'dd.mm.yyyy HH:MM:SS');
fprintf('Ended at %s\n\n', end_time);

end
