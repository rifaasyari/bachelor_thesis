function [s_estim] = ml_mimo(y, H, constellation, Possbl_s)
% ML_MIMO Implements the ML detection for MIMO systems
%   Input
%       y: Received signal
%       H: Channel gain matrix
%       constellation: Used transmit constellation, e.g. [-3, -1, 1, 3], 
%           [1, 1j, -1, -1j], ...
%
%   Output
%       s_estim: Estimation of transmitted signal
%
%   See also MB_LR_SIC

    N_t = cols(H);  % Number of transmit antennas of the system
    
    % from: https://de.mathworks.com/matlabcentral/fileexchange/40546-n-permute-k
    % Copyright (c) 2013, Adrian Etter
    % Copyright (c) 2006, Matt Fig
    % All rights reserved.
    % Possbl_s = npermutek(constellation, N_t)';  % Generate all possible sent 
                                                % signals for ML detection.
                                                % length(constellation)^N_t
                                                % permutations.
    % permutations = length(constellation)^N_t;   
    % addflops(flops_pow(N_t));
    
    Y = repmat(y, 1, cols(Possbl_s));
    ml_arg = Y - H*Possbl_s;
    ml_criterion = dot(ml_arg, ml_arg);
    
    % addflops(flops_mul(rows(H), N_t, permutations));
    % addflops(length(Y)); 
    % flops(flops + rows(ml_arg) * cols(ml_arg) + ... 
        % (rows(ml_arg) - 1) * cols(ml_arg));  
    % ADDFLOPS is sensitive for overflow. This prevents an overflow of FLOPS
  
    s_estim = Possbl_s(:, argmin(ml_criterion));
end
