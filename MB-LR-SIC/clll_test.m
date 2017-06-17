% CLLL Test

% % setup
% N_0 = 1;
% reps = 1000;
% trials = 20;
% N_t = 6;
% N_r = 6;
% E_b = 6:2:20;
% 
% vector_error_rate = zeros(1,20);
% for x = E_b
%     E_s = 8 * x;
%     constellation = create_MQAM(64, E_s);
%     
%     VER_mean = zeros(1,trials);
%     for i = 1:trials
%         H = 1/sqrt(2) * (randn(N_r,N_t) + 1j * randn(N_r,N_t));
%         
%         VER = 0;
%         for j = 1:reps
%             
%             s = constellation(randi([1,numel(constellation)], N_t, 1));
%             n = (sqrt(N_0) * randn(N_r,1)) .* exp(1j * 2 * pi * randn(N_r,1));
%             y = H * s + n;
% 
%             [H_LR, U] = clll(H, 0.99);
% 
%             % LR-aided SIC
%             [Q,R] = qr(H_LR);
%             z = ctranspose(Q) * y;
%             s_hat_LR = zeros(N_t,1);
%             for k = N_t:-1:1
%                 s_hat_LR(k) = z(k) / R(k,k);
%                 z = z - R(k, k+1:end) * s_hat_LR(k+1:end);
%             end
%             s_tilde = U*s_hat_LR;
%             % hard-limit componentwise to valid symbol vector
%             s_hat = zeros(N_t,1);
%             for k = 1:N_t
%                 S = repmat(s_tilde(k), 8, 8);
%                 ml_arg = abs(constellation-S).^2;
%                 [~, min_row] = min(ml_arg);
%                 [~, min_col] = min(min(ml_arg));
%                 s_hat(k) = constellation(min_row(min_col), min_col);
%             end
%             
%             if sum(s_hat ~= s)
%                 VER = VER + 1;
%             end
%         end
%         VER_mean(i) = VER / reps;
%         fprintf('E_b = %d, Trial = %d, VER = %.4f\n', x, i, VER_mean(i));
%     end
%     vector_error_rate(x) = mean(VER_mean);
% end
% 
% vector_error_rate = vector_error_rate(vector_error_rate > 1e-4);
% 
% figure
% semilogy(E_b, vector_error_rate, 'b--*', ...
%     'DisplayName', 'CLLL-SIC');
% xlabel('$E_b/N_0$ (dB)');
% ylabel('vector-error-rate');
% legend('show', 'Location', 'southwest');

N = 10^4;
n = 4;  % Number antennas
passed = 0;
for i = 1:N
    H = 1/sqrt(2) * (randn(n) + 1j*randn(n));
    [H_LR,T,P_c] = clll(H, 0.99);
    passed = passed + P_c;
end
fprintf('Probability that conditional test in line 106 is passed for n=%d: %.4f\n', ...
    n, passed / N);
