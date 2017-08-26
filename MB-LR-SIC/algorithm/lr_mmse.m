function [s_estim] = lr_mmse(y, H, signal_power, N, modulation, ...
    constellation, use_clll)

    % constellation = constellation(:);

    N_t = size(H, 2);  % Total number of transmit antennas is equal to 
                       % number of columns in H
    N_r = length(y);  
    
    if use_clll
        [H_LR, T] = clll(H, 0.99);  % CLLL Reduction 
        [~,R] = qr(H);
        [~,R_LR] = qr(H_LR);
        if norm(R-eye(N_t)) < norm(R_LR-eye(N_t))  % H is nearer to Orthogonal 
                                                   % than H_LR 
            H_LR = H;
            T = eye(N_r,N_t);
        end
    else
        H_LR = H;
        T = eye(size(H));
    end

    T_inv = inv(T); 

    R_zz = signal_power*inv(ctranspose(T)*T);

    % MMSE linear equalizer

    W = ctranspose(inv(H_LR*R_zz*ctranspose(H_LR)+N*eye(N_r))*...
                H_LR*R_zz);

    z_hat = W * y; % Output of MMSE detector. Detected signal component in 
                   % reduced lattice space

    % Shifting and Scaling operations 
    for n = 1:N_t
        z_hat(n) = shifting_scaling(z_hat(n), 1+1j, T_inv, n, ... 
            signal_power, modulation, false);
    end

    s_estim = T * z_hat;
        
    % Quantization: Return valid constellation symbol
    for i = 1:N_t
        S_L = repmat(s_estim(i), size(constellation));
        if ~any(S_L == constellation)
            errs = abs(S_L-constellation).^2;
            [~, argmin] = min(errs(:));
            s_estim(i) = constellation(argmin);
        end
    end
end