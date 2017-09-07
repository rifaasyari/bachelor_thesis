function [s_estim] = lr_sic(y, H, var_s, N, modulation, ...
    constellation, extended, use_clll)

    constellation = constellation(:)';  % Store constellation in a vector 
                                       % instead of a matrix


    N_t = size(H, 2);  % Total number of transmit antennas is equal to 
                       % number of columns in H
    N_r = size(H, 1);
    
    if nargin == 6
        extended = true;
    end
    if nargin <= 7
        use_clll = true;
    end
         
    if extended
        H = [H;sqrt(N)*eye(N_t)];
        y = [y;zeros(N_t,1)];
    end
    
    if use_clll
        [H_LR, T] = clll(H, 0.99);  % CLLL Reduction 
    else     
        % No CLLL
        H_LR = H;
        T = eye(N_t);
    end

    T_inv = inv(T); 

    % P_l = cnbo(H_LR,false);  % column-norm based ordering in H_LR
    P_l = pno(H_LR);

    H_LR_l = H_LR * P_l;  % Permute channel gain matrix for following SIC
    
    U = P_l'*T_inv;
    
    R_zz = var_s*P_l'*(T_inv*T_inv')*P_l;
    
    z_hat = zeros(N_t,1);
    y_l_n = y;

    for n = 1:N_t  % Loop over signal components

        % MMSE linear equalizer
        if ~extended
            R_yy = H_LR_l*R_zz*H_LR_l'+N*eye(N_r);
            R_zy = R_zz*H_LR_l';
            W_l_n = R_zy*inv(R_yy);
        else
            W_l_n = inv(H_LR_l'*H_LR_l)*H_LR_l';
        end

        z_l_n = W_l_n * y_l_n; % Output of MMSE detector. Detected 
                               % signal component in reduced lattice 
                               % space                                                

        if use_clll
            % Shifting and Scaling operations 
            z_hat(n) = shifting_scaling(z_l_n(1), 1+1j, U, n, ... 
                var_s, modulation);
        else
            Z = repmat(z_l_n(1), size(constellation));
            if ~any(Z == constellation)
                errs = abs(Z-constellation).^2;
                [~, argmin] = min(errs(:));
                z_hat(n) = constellation(argmin);
            else
                z_hat(n) = z_l_n(1);
            end
        end

        y_l_n = y_l_n - H_LR_l(:,1) * z_hat(n);
        
        H_LR_l = H_LR_l(:,2:end);  % Update channel gain matrix in LR 
                                   % domain by removing channel gains
                                   % for already estimated signal
                                   % components 
        R_zz = R_zz(2:end,2:end);
                                   
    end
    s_estim = T * P_l * z_hat;  % Estimation of transmitted signal for 
                                % the l-th branch
        
      
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
