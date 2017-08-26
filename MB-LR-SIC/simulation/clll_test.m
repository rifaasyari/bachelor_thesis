% Test results of CLLL

clear
clc

addpath('../algorithm');

% +++++++++++++ Probability that conditional test is passed ++++++++++++++

N = 10^3;
n = 4;  % Number antennas
passed = 0;
for i = 1:N
    H = 1/sqrt(2) * (randn(n) + 1j*randn(n));
    
    % My CLLL
    [H_LR,T,P_c] = clll(H, 0.99);
    passed = passed + P_c;
    

    assert(abs(det(T)) >= 1 - eps(1e4) && abs(det(T)) <= 1 + eps(1e4), 'T not uni-modular');
    assert(all(all(H_LR - H*T - eps(1e4) <= 0)), 'T no valid transformation');
    
    
end
fprintf('Probability that conditional test in line 106 is passed for n=%d: %.4f\n', ...
    n, passed / N);

% ++++++++++++++ Orthogonality of lattice reduced matrix ++++++++++++++++

trials = 1e3;
mat_size = 4:2:20;

orth_LR = zeros(max(mat_size),1);

for m = mat_size

    for i = 1:trials
        H = 1/sqrt(2) * (randn(m) + 1j*randn(m));
        H_LR = clll(H, 0.99);

        [~,R] = qr(H);
        [~,R_LR] = qr(H_LR);

        orth_R = norm(R - eye(m));
        orth_R_LR = norm(R_LR - eye(m));

        if orth_R_LR < orth_R
            orth_LR(m) = orth_LR(m) + 1;
        end
    end
    
end

orth_LR = orth_LR / trials;

% fprintf('Lattice reduced matrix is nearer to Orthogonal in %.2f%% ', ... 
%     orth_LR*100);
% fprintf('of %.1E random channel realization\n',trials);

disp(orth_LR);
