% BER simulations for the SD and FSD

addpath('../algorithm');

clear
clc

% +++++++++++++++++++ Begin configuration +++++++++++++++++++++++++++++

N_r = 8;  % Receive antennas
N_t = 8;  % Transmit antennas 

correlated_antennas = false;  % Consider spatial correlation in the channel

ml_dec = false && ~correlated_antennas;  % Use ML detector 
qam4 = true && ~correlated_antennas;  % Apply 4-QAM modulation
qam16 = true || correlated_antennas;  % Apply 16-QAM modulation

transmissions = 200;  % Transmissions per channel scenario
runs = 50;  % Channel scenarios per SNR

SNR_start = 0;
SNR_max = 20;  % SNR per bit defined as log2(P)^-1/N_0
SNR_step = 4;

% +++++++++++++++++++ End configuration ++++++++++++++++++++++++++++

E_s = 1/N_t;

if qam4 == true
    constellation_4 = create_MQAM(4,E_s);
    decode_4 = @decode_qpsk;
    constellation_4 = constellation_4(:);
    P_4 = numel(constellation_4);
end
if qam16 == true
    constellation_16 = create_16_QAM(E_s);
    decode_16 = @decode_16QAM;
    constellation_16 = constellation_16(:);
    P_16 = numel(constellation_16);
end

if ~correlated_antennas
    BER_sd_4 = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
    BER_fsd_4 = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
    BER_sd_16 = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
    BER_fsd_16 = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
    BER_ml = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
else
   BER_sd = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
   BER_sd_low = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
   BER_sd_medium = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
   BER_sd_high = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
   BER_fsd_low = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
   BER_fsd_medium = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
   BER_fsd_high = zeros(1, (SNR_max-SNR_start)/SNR_step+1);
   
   % Correlation matrices
   % Low correlation
   R_low = [1, 0.24-0.19j, 0.11+0.02j, 0.05+0.11j;
            0.24+0.19j, 1, 0.24-0.19j, 0.11+0.02j;
            0.11-0.02j, 0.24+0.19j, 1, 0.24-0.19j;
            0.05-0.11j, 0.11-0.02j, 0.24+0.19j, 1];
   % Medium correlation
   R_medium = [1, -0.50+0.05j, 0.21+0.11j, 0.01-0.11j;
              -0.50-0.05j, 1, -0.50+0.05j, 0.21+0.11j;
              0.21-0.11j, -0.50-0.05j, 1, -0.50+0.05j;
              0.01+0.11j, 0.21-0.11j, -0.50-0.05j, 1];
   % High correlation
   R_high = [1, 0.01+0.70j, -0.47-0.08j, 0.19-0.26j;
             0.01-0.70j, 1, 0.01+0.70j, -0.47-0.08j;
             -0.47+0.08j, 0.01-0.70j, 1, 0.01+0.70j;
             0.19+0.26j, -0.47+0.08j, 0.01-0.70j, 1];
end

opts = {'no', 'yes'};

% Print configuration
fprintf('++++++++++ Configuration ++++++++++\n');
fprintf('Transmit antennas: %d, Receive antennas: %d\n', N_t, N_r);
fprintf('ML detection: %s, 4-QAM: %s, 16-QAM: %s\n', ... 
    opts{ml_dec+1}, opts{qam4+1}, opts{qam16+1});
fprintf('Correlated antennas: %s\n', opts{correlated_antennas+1});
fprintf('Transmissions: %d, Simulation runs: %d\n', transmissions, runs);
fprintf('SNR: start = %d, step = %d, end = %d\n\n', SNR_start,SNR_step,SNR_max);

for SNR = SNR_start:SNR_step:SNR_max
    
    if qam4 == true
        N_0_4 = 1/log2(P_4)*10^(-SNR/10);
    end
    if qam16 == true
        N_0_16 = 1/log2(P_16)*10^(-SNR/10);
    end
    
    if ~correlated_antennas
        BERs_sd_4 = zeros(1,runs);
        BERs_fsd_4 = zeros(1,runs);
        BERs_sd_16 = zeros(1,runs);
        BERs_fsd_16 = zeros(1,runs);
        BERs_ml = zeros(1,runs);
    else
        BERs_sd = zeros(1,runs);
        BERs_sd_low = zeros(1,runs);
        BERs_sd_medium = zeros(1,runs);
        BERs_sd_high = zeros(1,runs);
        BERs_fsd_low = zeros(1,runs);
        BERs_fsd_medium = zeros(1,runs);
        BERs_fsd_high = zeros(1,runs);
    end
    
    for i = 1:runs  % channel realizations
        
        H = 1/sqrt(2) * (randn(N_r,N_t) + 1j*randn(N_r,N_t));
        
        if ~correlated_antennas
            BE_sd_4 = 0;
            BE_fsd_4 = 0;
            BE_sd_16 = 0;
            BE_fsd_16 = 0;
            BE_ml = 0;
        else
            BE_sd = 0;
            BE_sd_low = 0;
            BE_sd_medium = 0;
            BE_sd_high = 0;
            BE_fsd_low = 0;
            BE_fsd_medium = 0;
            BE_fsd_high = 0;
            
            % Correlated channel matrices
            H_low = sqrtm(R_low)*H*sqrtm(R_low);
            H_medium = sqrtm(R_medium)*H*sqrtm(R_medium);
            H_high = sqrtm(R_high)*H*sqrtm(R_high);
        end
        
        p = sym('P');
        l_p = ceil(sqrt(N_t)-1);  % Valid iff N_r == N_t
        l_1 = N_t - l_p;
        n_S = [ones(l_1,1); p*ones(l_p,1)];
        
        for j = 1:transmissions  % transmitted symbols
            
            if qam4 == true
                s_4 = constellation_4(randi([1,numel(constellation_4)], 1, N_t)); 
                n_4 = sqrt(N_0_4/2) * (randn(N_r,1) + 1j*randn(N_r,1));
                y_4 = H*s_4+n_4;
            end
            if qam16 == true
                s_16 = constellation_16(randi([1,numel(constellation_16)], 1, N_t)); 
                n_16 = sqrt(N_0_16/2) * (randn(N_r,1) + 1j*randn(N_r,1));
                y_16 = H*s_16+n_16;
                
                if correlated_antennas == true
                    y_low = H_low*s_16+n_16;
                    y_medium = H_medium*s_16+n_16;
                    y_high = H_high*s_16+n_16;
                end
            end
            
            if qam4 == true
                s_sd_4 = sd(y_4,H,N_t,constellation_4);
                s_fsd_4 = fsd(y_4,H,N_t,constellation_4);
                b_4 = decode_4(s_4,constellation_4);
                b_sd_4 = decode_4(s_sd_4,constellation_4);
                b_fsd_4 = decode_4(s_fsd_4,constellation_4);
                BE_sd_4 = BE_sd_4 + sum(b_sd_4 ~= b_4);
                BE_fsd_4 = BE_fsd_4 + sum(b_fsd_4 ~= b_4);
            end
            if qam16 == true
                s_sd_16 = sd(y_16,H,N_t,constellation_16);
                b_16 = decode_16(s_16,constellation_16);
                b_sd_16 = decode_16(s_sd_16,constellation_16);
                
                if ~correlated_antennas
                    s_fsd_16 = fsd(y_16,H,N_t,constellation_16);   
                    b_fsd_16 = decode_16(s_fsd_16,constellation_16);
                    BE_sd_16 = BE_sd_16 + sum(b_sd_16 ~= b_16);
                    BE_fsd_16 = BE_fsd_16 + sum(b_fsd_16 ~= b_16);
                else
                    s_sd_low = sd(y_low,H_low,N_t,constellation_16);
                    s_sd_medium = sd(y_medium,H_medium,N_t,constellation_16);
                    s_sd_high = sd(y_high,H_high,N_t,constellation_16);
                    s_fsd_low = fsd(y_low,H_low,N_t,constellation_16);
                    s_fsd_medium = fsd(y_medium,H_medium,N_t,constellation_16);
                    s_fsd_high = fsd(y_high,H_high,N_t,constellation_16);
                    
                    b_sd_low = decode_16(s_sd_low,constellation_16);
                    b_sd_medium = decode_16(s_sd_medium,constellation_16);
                    b_sd_high = decode_16(s_sd_high,constellation_16);
                    b_fsd_low = decode_16(s_fsd_low,constellation_16);
                    b_fsd_medium = decode_16(s_fsd_medium,constellation_16);
                    b_fsd_high = decode_16(s_fsd_high,constellation_16);
                    
                    BE_sd = BE_sd + sum(b_sd_16 ~= b_16);
                    BE_sd_low = BE_sd_low + sum(b_sd_low ~= b_16);
                    BE_sd_medium = BE_sd_medium + sum(b_sd_medium ~= b_16);
                    BE_sd_high = BE_sd_high + sum(b_sd_high ~= b_16);
                    BE_fsd_low = BE_fsd_low + sum(b_fsd_low ~= b_16);
                    BE_fsd_medium = BE_fsd_medium + sum(b_fsd_medium ~= b_16);
                    BE_fsd_high = BE_fsd_high + sum(b_fsd_high ~= b_16);
                end
            end
            if ml_dec == true
                [s_ml,err] = ml_mimo(y_4,H,constellation_4);
                b_ml = decode_4(s_ml,constellation_4);
                BE_ml = BE_ml + sum(b_ml ~= b_4);
            end 
        
        end
        
        if qam4 == true
            BERs_sd_4(i) = BE_sd_4 / (transmissions * N_t * ...
                log2(numel(constellation_4)));
            fprintf('%d) SNR=%d, 4-QAM, BER SD: %.4f\n', i, SNR, BERs_sd_4(i));
            BERs_fsd_4(i) = BE_fsd_4 / (transmissions * N_t * ...
                log2(numel(constellation_4)));
            fprintf('%d) SNR=%d, 4-QAM, BER FSD: %.4f\n', i, SNR, BERs_fsd_4(i));
        end
        
        if qam16 == true  && ~correlated_antennas
            BERs_sd_16(i) = BE_sd_16 / (transmissions * N_t * ...
                log2(numel(constellation_16)));
            fprintf('%d) SNR=%d, 16-QAM, BER SD: %.4f\n', i, SNR, BERs_sd_16(i));
            BERs_fsd_16(i) = BE_fsd_16 / (transmissions * N_t * ...
                log2(numel(constellation_16)));
            fprintf('%d) SNR=%d, 16-QAM, BER FSD: %.4f\n', i, SNR, BERs_fsd_16(i));
        elseif correlated_antennas
            BERs_sd(i) = BE_sd / (transmissions * N_t * ...
                log2(numel(constellation_16)));
            BERs_sd_low(i) = BE_sd_low / (transmissions * N_t * ...
                log2(numel(constellation_16)));
            BERs_sd_medium(i) = BE_sd_medium / (transmissions * N_t * ...
                log2(numel(constellation_16)));
            BERs_sd_high(i) = BE_sd_high / (transmissions * N_t * ...
                log2(numel(constellation_16)));
            BERs_fsd_low(i) = BE_fsd_low / (transmissions * N_t * ...
                log2(numel(constellation_16)));
            BERs_fsd_medium(i) = BE_fsd_medium / (transmissions * N_t * ...
                log2(numel(constellation_16)));
            BERs_fsd_high(i) = BE_fsd_high / (transmissions * N_t * ...
                log2(numel(constellation_16)));
            fprintf('SNR: %d, Pass: %d\n', SNR, i);
        end
        
        if ml_dec == true
            BERs_ml(i) = BE_ml / (transmissions * N_t * ... 
                log2(numel(constellation_4)));
            fprintf('%d) BER ML: %.4f\n', i, BERs_ml(i));
        end
        
    end
    
    x = (SNR-SNR_start)/SNR_step+1;
    if ~correlated_antennas
        BER_sd_4(x) = mean(BERs_sd_4);
        BER_fsd_4(x) = mean(BERs_fsd_4);
        BER_sd_16(x) = mean(BERs_sd_16);
        BER_fsd_16(x) = mean(BERs_fsd_16);
    else
        BER_sd(x) = mean(BERs_sd);
        BER_sd_low(x) = mean(BERs_sd_low);
        BER_sd_medium(x) = mean(BERs_sd_medium);
        BER_sd_high(x) = mean(BERs_sd_high);
        BER_fsd_low(x) = mean(BERs_fsd_low);
        BER_fsd_medium(x) = mean(BERs_fsd_medium);
        BER_fsd_high(x) = mean(BERs_fsd_high);
    end
    if ml_dec == true
        BER_ml((SNR-SNR_start)/SNR_step+1) = mean(BERs_ml);
    end
end

x_axis = SNR_start:SNR_step:SNR_max;

figure

if ~correlated_antennas
    semilogy(x_axis, BER_sd_4, 'k--^', ...
        'LineWidth', 1, ...
        'DisplayName', 'SD');
    hold on
    semilogy(x_axis, BER_fsd_4, 'k-s', ...
        'LineWidth', 1, ...
        'DisplayName', sprintf('FSD - n_S = ( %s)',sprintf('%s ',n_S)));

    if ml_dec == true
        hold on
        semilogy(x_axis, BER_ml, 'k-', ... 
                'LineWidth', 1, ...
                'DisplayName', 'ML');
    end

    legend('show', 'Location', 'northeast', 'AutoUpdate', 'off');

    if qam16 == true
        hold on
        semilogy(x_axis, BER_sd_16, 'k--^', ...
            'LineWidth', 1);
        hold on
        semilogy(x_axis, BER_fsd_16, 'k-s', ...
            'LineWidth', 1);
    end

    title(sprintf('$N_t=%d, N_r=%d$',N_t,N_r));
else
    semilogy(x_axis, BER_sd, 'k--', ...
        'LineWidth', 1, ...
        'DisplayName', 'SD, no correlation');
    hold on
    semilogy(x_axis, BER_sd_low, 'k--^', ...
        'LineWidth', 1, ...
        'DisplayName', 'SD');
    hold on
    semilogy(x_axis, BER_fsd_low, 'k-s', ...
        'LineWidth', 1, ...
        'DisplayName', sprintf('FSD - n_S = ( %s)',sprintf('%s ',n_S)));
    
    legend('show', 'Location', 'northeast', 'AutoUpdate', 'off');
    
    hold on
    semilogy(x_axis, BER_sd_medium, 'k--^', ...
        'LineWidth', 1);
    hold on
    semilogy(x_axis, BER_fsd_medium, 'k-s', ...
        'LineWidth', 1);
    hold on
    semilogy(x_axis, BER_sd_high, 'k--^', ...
        'LineWidth', 1);
    hold on
    semilogy(x_axis, BER_fsd_high, 'k-s', ...
        'LineWidth', 1);
    
    title(sprintf('$N_t=%d, N_r=%d$, 16-QAM',N_t,N_r));
end

xlabel('$E_b/N_0$ (dB)');
ylabel('BER');
ylim([0.9e-5, 1e0]);
xticks(SNR_start:SNR_step:SNR_max);
grid on

savefig(sprintf('BER_%s.fig', datestr(datetime('now'), 'dd-mm-yyyy_HH-MM-SS')));
close
