function [s_estim] = mb_lr_sic(y, H, var_s, N, modulation, ... 
    constellation, use_clll)
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
    
    constellation = constellation(:);  % Store constellation in a vector 
                                       % instead of a matrix

    if nargin == 6
        use_clll = true;
    end
    
    N_t = size(H, 2);  % Total number of transmit antennas is equal to 
                       % number of columns in H
    N_r = length(y);
    L = N_t;  % Number of permutations of channel gain matrix H
    
    H_orig = H;
    y_orig = y;
    
    H = [H;N*eye(N_t)];
    y = [y;zeros(N_t,1)];
    
    if use_clll
        [H_LR, T] = clll(H, 0.99);  % CLLL Reduction 
    else     
        % No CLLL
        H_LR = H;
        T = eye(N_t);
    end

    T_inv = inv(T); 

    % addflops(flops_inv(length(T)));
    
    S = zeros(N_t, L);
    for l = 1:L  % Multi-Branch loop
        if l == 1
            P_l = cnbo(H_LR);  % column-norm based ordering in H_LR % OK
        else
            P_l = psp(N_t, l, L);  % Pre_Stored Patterns % OK
        end
        
        H_LR_l = H_LR * P_l;  % Permute channel gain matrix for following 
                              % SIC
        % addflops(flops_mul(H_LR, P_l));
                
        z_hat = zeros(N_t,1);
        y_l_n = y;
        
        % T_l = T*P_l;  % Permuted transformation matrix
        
        % R_zz = var_s*inv(ctranspose(T_l)*T_l);
                                                        
        for n = 1:N_t  % Loop over signal components
                 
            % MMSE linear equalizer
%             W_l_n = ctranspose(inv(H_LR_l*R_zz*ctranspose(H_LR_l)+N*eye(N_r))*...
%                 H_LR_l*R_zz);
            W_l_n = inv(H_LR_l'*H_LR_l)*H_LR_l';
                        
            z_l_n = W_l_n * y_l_n; % Output of MMSE detector. Detected 
                                   % signal component in reduced lattice 
                                   % space    
                                                                                   
            % addflops(flops_mul(cols(W_l_n), rows(W_l_n), length(y_l_n)));

            % Shifting and Scaling operations 
            if use_clll
                z_hat(n) = shifting_scaling(z_l_n(1), 1+1j, P_l'*T_inv, n, ... 
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
            % addflops(2 * length(y_l_n));  
            
            H_LR_l = H_LR_l(:,2:end);  % Update channel gain matrix in LR 
                                       % domain by removing channel gains
                                       % for already estimated signal
                                       % components 
                                                                                          
        end
        
        s_l = T * P_l * z_hat;  % Estimation of transmitted signal for 
                                % the l-th branch   
      
        % Quantization: Return valid constellation symbol
        for i = 1:N_t
            S_L = repmat(s_l(i), size(constellation));
            if ~any(S_L == constellation)
                errs = abs(S_L-constellation).^2;
                [~, argmin] = min(errs(:));
                s_l(i) = constellation(argmin);
            end
        end
                                                       
        S(:,l) = s_l;
                                                                          
        % addflops(flops_mul(T, P_l) + flops_mul(N_t, N_t, 1));
        
    end
                                     
    Y = repmat(y_orig,1,L);  % Create matrix with L columns of received signal y 
                        % to quickly calulcate the best estimate of s
    ml_arg = Y - H_orig * S;
    [~, l_opt] = min(dot(ml_arg, ml_arg));  % ML decision
    
    s_estim = S(:,l_opt);
    
end 
