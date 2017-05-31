qpsk = [1*exp(0) 1*exp(pi/2 * 1j) 1*exp(pi * 1j) 1*exp(3*pi/2 * 1j)];
% qpsk modulation on the unit circle

antennas = 2:10;
reps=500;  % Compute the detection 'reps' times and then take average of flops
flops_mb_lr_sic = zeros(1, length(antennas));

for n = 1:length(antennas)
   N_t = antennas(n);
   N_r = N_t;  % same number of transmit and receive antennas
   flops_sum = 0;
   for i=1:reps
       H = randn(N_r, N_t) + 1j * randn(N_r, N_t);  % Compute random channel 
                                                    % gain matrix
       s = qpsk(randi([1,4], 1, N_t))';  % Generate signal vector
       y = H * s + (randn + 1j*randn);
       flops(0); mb_lr_sic(H, y, 1, 1); flops_sum = flops_sum + flops;
   end
   flops_mb_lr_sic(n) = flops_sum / reps;
end

plot(antennas, flops_mb_lr_sic);
