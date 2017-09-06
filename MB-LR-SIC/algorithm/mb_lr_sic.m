function [s_estim] = mb_lr_sic(y, H, var_s, N, modulation, ... 
    constellation, extended, use_clll)
% MB_LR_SIC Implements the "Mutli-Branch Lattice Reduction Successive
% Interference Cancellation Detection for Multiuser MIMO Systems" [1]
%
%   Input
%       y: Received signal 
%       H: Channel gain matrix
%       signal_power: Power of transmitted signals
%       N: Noise power
%       modulation (string): Used modulation: 'QPSK' or '16-QAM'
%       constellation (matrix):
%       use_clll (bool): If true apply CLLL reduction to H.
%           Otherwise this function performs *MB-SIC* detection
%   
%   Output
%       s_estim: Estimation of transmitted signal
%
%   [1]: L. Arvelo, R. Lamare, K. Zu, R. Sampaio-Neto, "Multi-Branch
%       Lattice Reduction Successive Interference Cancellation Detection 
%       for Multiuser MIMO Systems", 11th International Symposium on
%       Wireless Communication Systems, Barcelona, Spain, 2014, pp. 219 - 
%       223
%
%   See also CLLL, CNBO, PSP, MMSE_FILTER, SHIFTING_SCALING, ML_MIMO

    if nargin == 6
        extended = true;
    end
    if nargin <= 7
        use_clll = true;
    end
    
    N_t = size(H, 2);  % Total number of transmit antennas is equal to 
                       % number of columns in H
    N_r = size(H, 1);
    L = N_t + 1;  % Number of permutations of channel gain matrix H
    
    H_orig = H;
    y_orig = y;
    
    if extended
        H = [H;N*eye(N_t)];
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
    
    S = zeros(N_t, L);
    for l = 1:L  % Multi-Branch loop
        if l == 1
            P_l = cnbo(H_LR,false);  % column-norm based ordering in H_LR
        elseif l == 2
            P_l = pno(H_LR);
        else
            P_l = psp(N_t, l, L);  % Pre_Stored Patterns
        end
        
        H_LR_l = H_LR * P_l;  % Permute channel gain matrix for following 
                              % SIC
         
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
                                                                                   
            % Shifting and Scaling operations 
            if use_clll
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

            y_l_n = y_l_n - H_LR_l(:,1) * z_hat(n);  % Subtract estimated 
                                                     % component from
                                                     % received signal to
                                                     % further reduce
                                                     % interference for
                                                     % detection in next
                                                     % iteration
            
            H_LR_l = H_LR_l(:,2:end);  % Update channel gain matrix in LR 
                                       % domain by removing channel gains
                                       % for already estimated signal
                                       % components 
            R_zz = R_zz(2:end,2:end);
        end
        
        S(:,l) = T * P_l * z_hat;  % Estimation of transmitted signal for 
                                % the l-th branch
      
        % Quantization: Return valid constellation symbol
        for i = 1:N_t
            S_L = repmat(S(i,l), size(constellation));
            if ~any(S_L == constellation)
                errs = abs(S_L-constellation).^2;
                [~, argmin] = min(errs(:));
                S(i,l) = constellation(argmin);
            end
        end
                                                                                                                                         
    end
                                     
    Y = repmat(y_orig,1,L);  % Create matrix with L columns of received signal 
                             % y to quickly calulcate the best estimate of s
    ml_arg = Y - H_orig * S;
    [~, l_opt] = min(dot(ml_arg, ml_arg));  % ML decision
    
    s_estim = S(:,l_opt);
    
end 
