qpsk = [1*exp(0) 1*exp(pi/2 * 1j) 1*exp(pi * 1j) 1*exp(3*pi/2 * 1j)];
% qpsk modulation on the unit circle

antennas = 2:10;
flops_mb_lr_sic = zeros(1, length(antennas));
flops_ml = zeros(1, length(antennas));

for n = 1:length(antennas)
   N_t = antennas(n);
   N_r = N_t;  % same number of transmit and receive antennas
   
   H = randn(N_r, N_t) + 1j * randn(N_r, N_t);  % Compute random channel 
                                               % gain matrix
   s = qpsk(randi([1,4], 1, N_t))';  % Generate signal vector
   y = H * s + (randn(N_r,1) + 1j*randn(N_r,1));
   flops(0); mb_lr_sic(y, H, 1, 'QPSK'); flops_mb_lr_sic(n) = flops;
   flops(0); ml_mimo(y, H, qpsk); flops_ml(n) = flops;
   
end

% plot
figure
semilogy(antennas, flops_mb_lr_sic, 'r-s', ... 
    'LineWidth', 2, ...
    'DisplayName', 'MB-LR-SIC');
hold on
semilogy(antennas, flops_ml, 'k-+', ... 
    'LineWidth', 2, ...
    'DisplayName', 'ML-QPSK');
title('Computational Complexity for MIMO System, L=$N_t$');
xlabel('Number of Antennas ($N_t$=$N_r$)');
ylabel('Number of Flops');
legend('show', 'Location', 'southeast');
