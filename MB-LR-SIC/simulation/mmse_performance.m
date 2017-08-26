% Performance of the MMMSE linear filter with and without lattice reduction
% (LR). No successive interference cancellation (SIC) is considered

clear 
clc

addpath('../algorithm');

% ++++++++++++++++++++++ Begin configuration +++++++++++++++++++++++++++

trials = 1e4;
N_t = 4;
N_r = N_t;
modulation = 'QPSK';
N = 0.08;  % N = 1 -> SNR = 9dB, N = 0.08 -> SNR = 20 dB
SS = true;  % Perform shifting and scaling operations for LR signal
slicing = true;  % Slice not-LR-signal to QAM table
extended = true;  % Use extended channel matrix and receive vector

% +++++++++++++++++++++++ End configuration ++++++++++++++++++++++++++++

if strcmp(modulation,'QPSK')
    E_s = 2;
    constellation = create_MQAM(4,E_s);
elseif strcmp(modulation,'16-QAM')
    E_s = 2;
    constellation = create_MQAM(16,E_s);
else
    error('Invalid modulation: %s',modulation);
end

disp(constellation);

SNR = 10*log10(N_t*E_s/N);

err_s = 0;
err_z = 0;
err_s_from_z = 0;

err_s_1st_addend = 0;
n_reduc_s = 0;

err_z_1st_addend = 0;
n_reduc_z = 0;

err_s_from_z_1st_addend = 0;

for i = 1:trials
    H = 1/sqrt(2) * (randn(N_r,N_t) + 1j*randn(N_r,N_t));
    
    s = constellation(randi([1,numel(constellation)], 1, N_t))';
    n = sqrt(N/2) * (randn(N_r,1) + 1j*randn(N_r,1));

    y = H*s+n;
    
    if extended == true
        H = [H; N*eye(N_t)];
        y = [y; zeros(N_t,1)];
        n = [n; zeros(N_t,1)]; 
    end
    [H_LR,T] = clll(H,0.99);
    T_inv = inv(T);
    
    z = T_inv*s;
       
    if ~extended
        W = E_s*H'*inv(E_s*(H*H')+N*eye(N_r));
        R_zz = E_s*(T_inv*T_inv');
        R_yy = H_LR*R_zz*H_LR'+N*eye(N_r);
        R_zy = R_zz*H_LR';
        W_LR = R_zy*inv(R_yy);
    else
        W = inv(H'*H)*H';
        W_LR = inv(H_LR'*H_LR)*H_LR';
    end
    
    s_tilde = W*y;
    if slicing == true
        for j = 1:N_t
            S = repmat(s_tilde(j),size(constellation));
            if ~any(S == constellation)
                errs = abs(S-constellation).^2;
                [~, argmin] = min(errs(:));
                s_tilde(j) = constellation(argmin);
            end
        end
    end
    z_tilde = W_LR*y;
    if SS == true
        if strcmp(modulation, '16-QAM')
            d = sqrt(2/5 * E_s);
        else
            d = sqrt(2 * E_s);
        end
        alpha = 2/d;
        z_tilde = z_tilde * alpha;
        z_tilde = 2*round((z_tilde-(1+1j)*sum(T_inv,2))/2) + ... 
            (1+1j)*sum(T_inv,2);
        z_tilde = z_tilde / alpha;
    end
    s_from_z = T*z_tilde;
    if slicing == true
        for j = 1:N_t
            S = repmat(s_from_z(j),size(constellation));
            if ~any(S == constellation)
                errs = abs(S-constellation).^2;
                [~, argmin] = min(errs(:));
                s_from_z(j) = constellation(argmin);
            end
        end
    end
    
    est_s_1st_addend = W*(H*s);
    if slicing == true
        for j = 1:N_t
            S = repmat(est_s_1st_addend(j),size(constellation));
            if ~any(S == constellation)
                errs = abs(S-constellation).^2;
                [~, argmin] = min(errs(:));
                est_s_1st_addend(j) = constellation(argmin);
            end
        end
    end
    est_z_1st_addend = W_LR*(H_LR*z);
    if SS == true
        if strcmp(modulation, '16-QAM')
            d = sqrt(2/5 * E_s);
        else
            a = sqrt(2 * E_s);
        end
        alpha = 2/d;
        est_z_1st_addend = est_z_1st_addend * alpha;
        est_z_1st_addend = 2*round((est_z_1st_addend-(1+1j)*sum(T_inv,2))/2)+...
            (1+1j)*sum(T_inv,2);
        est_z_1st_addend = est_z_1st_addend / alpha;
    end
    
    s_from_z_1st_addend = T*est_z_1st_addend;
    if slicing == true
        for j = 1:N_t
            S = repmat(s_from_z_1st_addend(j),size(constellation));
            if ~any(S == constellation)
                errs = abs(S-constellation).^2;
                [~, argmin] = min(errs(:));
                s_from_z_1st_addend(j) = constellation(argmin);
            end
        end
    end
    
    est_s_n = W*n;
    est_z_n = W_LR*n;
    
    err_s = err_s + norm(s-s_tilde);
    err_z = err_z + norm(z-z_tilde);
    err_s_from_z = err_s_from_z + norm(s-s_from_z);
    
    err_s_1st_addend = err_s_1st_addend + norm(est_s_1st_addend-s);
    err_z_1st_addend = err_z_1st_addend + norm(est_z_1st_addend-z);
    err_s_from_z_1st_addend = err_s_from_z_1st_addend + ...
        norm(s_from_z_1st_addend-s);
    
    n_reduc_s = n_reduc_s + (norm(est_s_n)/norm(n));
    n_reduc_z = n_reduc_z + (norm(est_z_n)/norm(n));
end

err_s = err_s / trials;
err_z = err_z / trials;
err_s_from_z = err_s_from_z / trials;

err_s_1st_addend = err_s_1st_addend / trials;
err_z_1st_addend = err_z_1st_addend / trials;
err_s_from_z_1st_addend = err_s_from_z_1st_addend / trials;

n_reduc_s = n_reduc_s / trials;
n_reduc_z = n_reduc_z / trials;

fprintf('SNR = %d dB\n', round(SNR));

if extended == true
    fprintf('Use extended channel matrix and receive vector\n\n');
else
    fprintf('\n');
end

if ~slicing
    fprintf('Average error in s (w/o QAM slicing) after %.1E trials: %.4f\n',...
        trials,err_s);
else
    fprintf('Average error in s (with QAM slicing) after %.1E trials: %.4f\n',...
        trials,err_s);
end
if ~SS
    fprintf('Average error in z (w/o SS) after %.1E trials: %.4f\n',... 
        trials,err_z);
else
    fprintf('Average error in z (with SS) after %.1E trials: %.4f\n',... 
        trials,err_z);
end
if ~slicing && ~SS
    fprintf('Average error in s estimated from z (w/o SS and QAM slicing) after ');
    fprintf('%.1E trials: %.4f\n\n', trials,err_s_from_z);
elseif ~slicing && SS
    fprintf('Average error in s estimated from z (with SS and w/o QAM slicing) after ');
    fprintf('%.1E trials: %.4f\n\n', trials,err_s_from_z);
elseif slicing && ~SS
    fprintf('Average error in s estimated from z (w/o SS and with QAM slicing) after ');
    fprintf('%.1E trials: %.4f\n\n', trials,err_s_from_z);
else
    fprintf('Average error in s estimated from z (with SS and QAM slicing) after ');
    fprintf('%.1E trials: %.4f\n\n', trials,err_s_from_z);
end

if ~slicing
    fprintf('Average error in 1st addend of s (w/o QAM slicing) after %.1E trials: %.4f\n', ...
        trials, err_s_1st_addend);
else
    fprintf('Average error in 1st addend of s (with QAM slicing) after %.1E trials: %.4f\n', ...
        trials, err_s_1st_addend);
end
if ~SS
    fprintf('Average error in 1st addend of z (w/o SS) after %.1E trials: %.4f\n', ...
        trials, err_z_1st_addend);
else
    fprintf('Average error in 1st addend of z (with SS) after %.1E trials: %.4f\n', ...
        trials, err_z_1st_addend);
end
if ~slicing && ~SS
    fprintf('Average error in 1st addend of s estimated from z ');
    fprintf('(w/o SS and QAM slicing) after %.1E trials: %.4f\n\n', ... 
        trials,err_s_from_z_1st_addend);
elseif ~slicing && SS
    fprintf('Average error in 1st addend of s estimated from z ');
    fprintf('(with SS and w/o QAM slicing) after %.1E trials: %.4f\n\n', ... 
        trials,err_s_from_z_1st_addend);
elseif slicing && ~SS
    fprintf('Average error in 1st addend of s estimated from z ');
    fprintf('(w/o SS and with QAM slicing) after %.1E trials: %.4f\n\n', ... 
        trials,err_s_from_z_1st_addend);
else
    fprintf('Average error in 1st addend of s estimated from z ');
    fprintf('(with SS and QAM slicing) after %.1E trials: %.4f\n\n', ... 
        trials,err_s_from_z_1st_addend);
end

fprintf('Ratio ||W*n||/||n|| after %.1E trials: %.4f\n',trials,n_reduc_s);
fprintf('Ratio ||W_LR*n||/||n|| after %.1E trials: %.4f\n\n',trials,n_reduc_z);
