% ++++++++++++++++++++++ Complexity of SD ++++++++++++++++++++++++++++

addpath('../algorithm');

clear
clc

N_r = 4;
N_t = 4;
SNR = 0;

% 4-QAM
E_s_4 = 2;
N_0_4 = E_s_4 * 10^(-SNR/10);
constellation_4 = create_MQAM(4,E_s_4);
decode_4 = @decode_qpsk;
constellation_4 = constellation_4(:);

% 16-QAM
E_s_16 = 10;
N_0_16 = E_s_16 * 10^(-SNR/10);
constellation_16 = create_16_QAM(E_s_16);
decode_16 = @decode_16QAM;
constellation_16 = constellation_16(:);

transmissions = 200;

paths_4 = zeros(1,transmissions);
all_4 = zeros(1,transmissions);
paths_16 = zeros(1,transmissions);
all_16 = zeros(1,transmissions);

for i = 1:transmissions
    H = 1/sqrt(2)*(randn(N_r,N_t) + 1j*randn(N_r,N_t));
      
    s_16 = constellation_4(randi([1,numel(constellation_4)], 1, N_t)); 
    n_16 = sqrt(N_0_4/2) * (randn(N_r,1) + 1j*randn(N_r,1));
    y_4 = H*s_16+n_16;
    
    s_16 = constellation_16(randi([1,numel(constellation_16)], 1, N_t)); 
    n_16 = sqrt(N_0_16/2) * (randn(N_r,1) + 1j*randn(N_r,1));
    y_16 = H*s_16+n_16;
    
    [s_sd_4,R_4,~,~,full_path_4, visited_4] = sd(y_4,H,N_t,constellation_4);
    [s_sd_16,R_16,~,~,ful_path_16, visited_16] = sd(y_16,H,N_t,constellation_16);
    
    paths_4(i) = full_path_4;
    all_4(i) = visited_4;
    paths_16(i) = ful_path_16;
    all_16(i) = visited_16;
    
end

figure
bar(1:transmissions, paths_16, 'g', ...
    'DisplayName', '16-QAM');
hold on
bar(1:transmissions, paths_4, 'r', ...
    'DisplayName', '4-QAM');

title(sprintf('$N_r=N_t=%d, SNR = %d dB$', N_t, SNR));
xlabel('Transmissions');
xlim([0,transmissions]);
ylabel('Fully Searched Paths');
legend('show', 'location', 'northeast');

% ++++++++++++++ Time comparaison between SD and FSD +++++++++++++++++++++

N_t = 4;
N_r = N_t;

transmissions = 1e3;

E_s = 1/N_t;
% 16-QAM
constellation_16 = create_16_QAM(E_s);
constellation_16 = constellation_16(:);
P_16 = numel(constellation_16);

% 64-QAM
constellation_64 = create_MQAM(64,E_s);
constellation_64 = constellation_64(:);
P_64 = numel(constellation_64);

p = sym('P');
l_p = ceil(sqrt(N_t)-1);  % Valid iff N_r == N_t
l_1 = N_t - l_p;
n_S = [ones(l_1,1); p*ones(l_p,1)];

SNR_start = 0;
SNR_max = 30;
SNR_step = 5;

time_sd_16= zeros(1,(SNR_max-SNR_start)/SNR_step+1);
time_fsd_16 = zeros(1,(SNR_max-SNR_start)/SNR_step+1);
time_sd_64 = zeros(1,(SNR_max-SNR_start)/SNR_step+1);
time_fsd_64 = zeros(1,(SNR_max-SNR_start)/SNR_step+1);

i = 1;

for SNR = SNR_start:SNR_step:SNR_max
    N_0_16 = 1/log2(P_16)*10^(-SNR/10);
    N_0_64 = 1/log2(P_64)*10^(-SNR/10);
    for j = 1:transmissions
        H = 1/sqrt(2) * (randn(N_r,N_t) + 1j*randn(N_r,N_t));
        
        s_16 = constellation_16(randi([1,P_16], 1, N_t)); 
        n_16 = sqrt(N_0_16/2) * (randn(N_r,1) + 1j*randn(N_r,1));
        y_16 = H*s_16+n_16;
        
        s_64 = constellation_64(randi([1,P_64], 1, N_t));
        n_64 = sqrt(N_0_64/2) * (randn(N_r,1) + 1j*randn(N_r,1));
        y_64 = H*s_64+n_64;
        
        % Time for 16-QAM
        tic
            s_sd_16 = sd(y_16,H,N_t,constellation_16);
        t = toc;
        time_sd_16(i) = time_sd_16(i) + t;
        
        tic
            s_fsd_16 = fsd(y_16,H,N_t,constellation_16);
        t = toc;
        time_fsd_16(i) = time_fsd_16(i) + t;
        
        % Time for 64-QAM
        tic
            s_sd_64 = sd(y_64,H,N_t,constellation_64);
        t = toc;
        time_sd_64(i) = time_sd_64(i) + t;
        
        tic
            s_fsd_64 = fsd(y_64,H,N_t,constellation_64);
        t = toc;
        time_fsd_64(i) = time_fsd_64(i) + t;
        
    end
    i = i + 1;
end

time_sd_16 = time_sd_16 / transmissions;
time_fsd_16 = time_fsd_16 / transmissions;
time_sd_64 = time_sd_64 / transmissions;
time_fsd_64 = time_fsd_64 / transmissions;

figure
plot(SNR_start:SNR_step:SNR_max,time_sd_16, 'k-^', ...
    'LineWidth', 1, ...
    'DisplayName', 'SD');
hold on
plot(SNR_start:SNR_step:SNR_max,time_fsd_16, 'k-s', ...
    'LineWidth', 1, ...
    'DisplayName', sprintf('FSD - n_S = ( %s)',sprintf('%s ',n_S)));

legend('show', 'location', 'northeast', ...
    'AutoUpdate', 'off');

hold on
plot(SNR_start:SNR_step:SNR_max,time_sd_64, 'k-^', ...
    'LineWidth', 1);
hold on
plot(SNR_start:SNR_step:SNR_max,time_fsd_64, 'k-s', ...
    'LineWidth', 1);

title(sprintf('$N_t=N_r=%d$',N_t));
xlabel('$E_b/N_0$ (dB)');
ylabel('Time / s');
xticks(SNR_start:SNR_step:SNR_max)
ylim([0 0.06]);

grid on;
