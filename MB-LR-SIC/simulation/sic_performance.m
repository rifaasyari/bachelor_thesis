% Performance of the MB-SIC-MMMSE with and without lattice reduction (LR). 

clear 
clc

addpath('../algorithm');

% ++++++++++++++++++++++ Begin configuration +++++++++++++++++++++++++++

trials = 1e4;
N_t = 4;
N_r = N_t;
L = N_t+1;  % 1 <= L <= N_t
modulation = '16-QAM';
E_s = 1;
N = 1;  % N = 1 -> SNR = 9dB, N = 0.08 -> SNR = 20dB
extended = true;
QR_decomp = false && extended;

% +++++++++++++++++++++++ End configuration ++++++++++++++++++++++++++++

if strcmp(modulation,'QPSK')
    constellation = create_MQAM(4,E_s);
elseif strcmp(modulation,'16-QAM')
    constellation = create_MQAM(16,E_s);
else
    error('Invalid modulation: %s\n', modulation);
end

SNR = 10*log10(N_t*E_s/N);

err_s = 0;
err_z = 0;
err_s_from_z = 0;

for i = 1:trials
    H = 1/sqrt(2) * (randn(N_r,N_t) + 1j*randn(N_r,N_t));
    
    s = constellation(randi([1,numel(constellation)], 1, N_t))';
    n = sqrt(N/2) * (randn(N_r,1) + 1j*randn(N_r,1));
    
    y = H*s+n;
    
    H_orig = H;
    y_orig = y;
    
    if extended == true
        H = [H; N*eye(N_t)];
        y = [y; zeros(N_t,1)];
        % n = [n; zeros(N_t,1)]; 
    end
    
    [H_LR,T] = clll(H,0.99);
    [H_LR_orig,T_orig] = clll(H_orig,0.99);
    T_inv = inv(T);
    
    Z = repmat(inv(T_orig)*s,1,L);
    
    s_hat = zeros(N_t,L);
    z_hat = zeros(N_t,L);
    s_from_z = zeros(N_t,L);
        
    for l = 1:L
    
        if l == 1
            P = cnbo(H);
            P_LR = cnbo(H_LR);
        elseif l == 2
            P = sqrd(H);
            P_LR = sqrd(H_LR);
        else
            P = psp(N_t,l,L);
            P_LR = P;
        end
            
        H_P = H*P;
        H_LR_P = H_LR*P_LR;
        
        if QR_decomp
            [Q,R] = qr(H_P);
            [Q_LR,R_LR] = qr(H_LR_P);
        end

        U = P_LR'*T_inv;

        z = P_LR'*T_inv*s;
        
        R_zz = E_s*P_LR'*(T_inv*T_inv')*P_LR;
        
        if QR_decomp
            y1 = Q'*y;
            y2 = Q_LR'*y;
            dec_order = N_t:-1:1;
        else
            y1 = y;
            y2 = y;
            dec_order = 1:N_t;
        end

        for k = dec_order
            
            if ~extended
                W = E_s*H_P'*inv(E_s*(H_P*H_P')+N*eye(N_r));
                R_yy = H_LR_P*R_zz*H_LR_P'+N*eye(N_r);
                R_zy = R_zz*H_LR_P';
                W_LR = R_zy*inv(R_yy);
            else
                W = inv(H_P'*H_P)*H_P';
                W_LR = inv(H_LR_P'*H_LR_P)*H_LR_P';
            end

            if QR_decomp
                s_tilde = (y1(k)-R(k,k+1:end)*s_hat(k+1:end,l))/R(k,k);
            else
                s_tilde = W*y1;
            end
            S = repmat(s_tilde(1),size(constellation));
            % QAM slicing
            if ~any(S == constellation)
                errs = abs(S-constellation).^2;
                [~, argmin] = min(errs(:));
                s_hat(k,l) = constellation(argmin);
            end             

            if QR_decomp
                z_tilde = (y2(k)-R_LR(k,k+1:end)*z_hat(k+1:end,l))/R_LR(k,k);
            else
                z_tilde = W_LR*y2;
            end
            % Shifting and scaling
            z_hat(k,l) = shifting_scaling(z_tilde(1), 1+1j, U, k, ... 
                    E_s, modulation);

            if ~QR_decomp
                y1 = y1 - H_P(:,1)*s_hat(k,l);
                H_P = H_P(:,2:end);

                y2 = y2 - H_LR_P(:,1)*z_hat(k,l);
                H_LR_P = H_LR_P(:,2:end);
                R_zz = R_zz(2:end,2:end);
            end
        end
        
        s_hat(:,l) = P*s_hat(:,l);
        s_from_z(:,l) = T*P_LR*z_hat(:,l);
        z_hat(:,l) = P_LR*z_hat(:,l);
        
    end
    
    S = repmat(s,1,L);
    Y_orig = repmat(y_orig,1,L);
    Y = repmat(y(1:N_t),1,L);
    
    err_s = err_s + min(sqrt(dot(Y_orig-H_orig*s_hat,Y_orig-H_orig*s_hat)));
    err_z = err_z + min(sqrt(dot(Y-H_LR(1:N_t,:)*z_hat,Y-H_LR(1:N_t,:)*z_hat)));
    err_s_from_z = err_s_from_z + ... 
        min(sqrt(dot(Y_orig-H_orig*s_from_z,Y_orig-H_orig*s_from_z)));
end

err_s = err_s / trials;
err_z = err_z / trials;
err_s_from_z = err_s_from_z / trials;

fprintf('SNR = %d dB\n',round(SNR));

if extended
    fprintf('Use extendend channel matrix and receive vector\n');
end

if QR_decomp
    fprintf('Use QR decomposition for SIC\n');
end

fprintf('Average error in s after %.1E trials: %.4f\n',trials,err_s);

fprintf('Average error in z after %.1E trials: %.4f\n',trials,err_z);

fprintf('Average error in s estimated from z after %.1E trials: %.4f\n\n',...
    trials,err_s_from_z);
