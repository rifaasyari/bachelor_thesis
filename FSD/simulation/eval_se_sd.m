% Average number of nodes visited per parent node per level in SE-SD
% nodes denotes "child nodes visited per parent node"

addpath('../algorithm');

clear
clc

% The paper evaluated the average number of visited nodes per level for a
%   4x4 system with 16-QAM modulation
N_r = 4;
N_t = 4;
E_s = 1/N_t;

constellation = create_16_QAM(E_s);
constellation = constellation(:);
P = numel(constellation);

trials = 1e3;

SNR = 0:10:20;

average_n = zeros(N_t, length(SNR));
col = 1;

for snr = SNR
    N_0 = 1/log2(P)*10^(-snr/10);
    
    for i = 1:trials
        H = 1/sqrt(2) * (randn(N_r,N_t) + 1j*randn(N_r,N_t));
        % H = fsd_ordering(H,16,[1,1,1,16]);
        s = constellation(randi([1,numel(constellation)], 1, N_t)); 
        n = sqrt(N_0/2) * (randn(N_r,1) + 1j*randn(N_r,1));
        y = H*s+n;
        
        [s_sd,R,E_n] = sd(y,H,N_t,constellation);
        average_n(:,col) = average_n(:,col) + E_n;
    end
    
    col = col + 1;
end

average_n = average_n / trials;

disp(average_n);
