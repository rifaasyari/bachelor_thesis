function [s_estim] = mb_lr_sic(H, y, signal_power, noise_power)
% MB_LR_SIC Implements the "Mutli-Branch Lattice Reduction Successive
% Interference Cancellation Detection for Multiuser MIMO Systems" [1]
%
%   Input
%       H: Channel gain matrix
%       y: Received signal 
%       signal_power: Power of transmitted signals
%       noise_power: Power of AWGN noise over channel
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
%   See also CLLL, CNBO, PSP, MMSE_FILTER, SHIFTING_SCALING

    N_t = size(H, 2);  % Total number of transmit antennas is equal to 
                       % number of columns in H
    L = N_t;  % Number of permutations of channel gain matrix H
    [H_LR, T] = clll(H, 1);  % CLLL Reduction
    T_inv = inv(T);
    
    S = zeros(L);
    for l = 1:L  % Multi-Branch loop
        if l == 1
            P_l = cnbo(H_LR);  % column-norm based ordering in H_LR
        else
            P_l = psp(N_t, l, L);  % Pre_Stored Patterns
        end
        
        H_LR_l = H_LR * P_l;  % Permutate channel gain matrix for following 
                              % SIC
        
        z_estim = zeros(N_t);
        y_l_n = y;
        for n = 1:N_t  % Loop over signal components
            W_l_n = mmse_filter(H_LR_l, signal_power, noise_power); % MMSE 
            % linear equalizer
            z_l_n = ctranspose(W_l_n) * y_l_n;  % Output of MMSE detector. 
                                                % Detected signal component
                                                % in reduced lattice space
            z_estim(n) = shifting_scaling(z_l_n, n, 1+1j, T_inv);  % Shifting 
            % and Scaling operations 
            y_l_n = y_l_n - H_LR_l(:,1) * z_estim(n);  % Subtract estimated 
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
        end
        S(:,l) = T * P_l * z_estim;  % Estimation of transmitted signal for 
                                     % the l-th branch
    end
    Y = repmat(y,1,L);  % Create matrix with L columns of received signal y 
                        % to quickly calulcate the best estimate of s
    [err, l_opt] = min(sum((Y - H * S)^2, 1));  % ML decision
    
    s_estim = S(:, l_opt);
end 
