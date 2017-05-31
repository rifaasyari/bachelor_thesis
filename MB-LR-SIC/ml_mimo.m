function [s_estim] = ml_mimo(y, H, constellation)
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
    Possbl_s = npermutek(constellation, N_t)';  % Generate all possible sent 
                                                % signals for ML detection
                                               
    addflops(flops_pow(K));
    
    Y = repmat(y, 1, cols(Possbl_s));
    ml_criterion = dot(Y - H*Possbl_s);
    
    addflops();
    
    s_estim = Possbl_s(:, argmin(ml_criterion));
end
