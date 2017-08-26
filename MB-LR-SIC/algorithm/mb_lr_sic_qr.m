function [s_estim] = mb_lr_sic_qr(y, H, var_s, N, modulation, ... 
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
    
    if use_clll
        [H_LR, T] = clll(H, 0.99);  % CLLL Reduction 
        [~,R] = qr(H);
        [~,R_LR] = qr(H_LR);
        if norm(R - eye(N_t)) < norm(R_LR - eye(N_t))
            H_LR = H;
            T = eye(N_r,N_t);
        end
    else     
        % No CLLL
        H_LR = H;
        T = eye(size(H));
    end

    T_inv = inv(T); 

    % addflops(flops_inv(length(T)));
    
    S = zeros(N_t, L);
    for l = 1:L  % Multi-Branch loop
        if l == 1
            P_l = cnbo(H_LR);  % column-norm based ordering in H_LR
        else
            P_l = psp(N_t, l, L);  % Pre_Stored Patterns
        end
        
        H_LR_l = H_LR * P_l;  % Permute channel gain matrix for following 
                              % SIC

        % addflops(flops_mul(H_LR, P_l));
        
        [Q, R] = qr(H_LR_l);
                
        z_hat = zeros(N_t,1);
        y_l_n = ctranspose(Q) * y;
        
        T_l = T*P_l;  % Permuted transformation matrix
        
        R_zz = inv(ctranspose(T_l)*T_l);
                                        
        for n = 1:N_t  % Loop over signal components
                        
            % MMSE linear equalizer
            W_l_n = ctranspose(inv(R*R_zz*ctranspose(R)+N*eye(N_r))*R*R_zz);
            
            z_l_n = W_l_n * y_l_n; % Output of MMSE detector. Detected 
                                   % signal component in reduced lattice 
                                   % space                                                
                                                
            % addflops(flops_mul(cols(W_l_n), rows(W_l_n), length(y_l_n)));

            % Shifting and Scaling operations 
            z_hat(n) = shifting_scaling(z_l_n(1), 1+1j, T_inv, P_l, n, ... 
                var_s, modulation, true);

            y_l_n = y_l_n - R(:,1) * z_hat(n);

            % addflops(2 * length(y_l_n));  
                        
            R = R(:,2:end);
            
            R_zz = R_zz(2:end,2:end);
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
    
    % [Q, R] = qr(H);
                                     
    Y = repmat(y,1,L);  % Create matrix with L columns of received signal y 
                        % to quickly calulcate the best estimate of s
    ml_arg = Y - H * S;
    [~, l_opt] = min(dot(ml_arg, ml_arg));  % ML decision
    
    s_estim = S(:,l_opt);
    
end 
