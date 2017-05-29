function [s_estim] = ml_mimo(y, H, modulation)
% ML_MIMO Implements the ML detection for MIMO systems
%   Input
%       y: Received signal
%       H: Channel gain matrix
%       modulation: Used transmit modulation, e.g. M-QAM, QPSK, ...
%
%   Output
%       s_estim: Estimation of transmitted signal
%
%   See also MB_LR_SIC

    N_t = size(H, 2);  % Number of transmit antennas of the system
    
    % from: https://stackoverflow.com/questions/18591440/how-to-find-all-permutations-with-repetition-in-matlab
    C = cell(N_t, 1);
    [C{:}] = ngrid(modulation);
    possible_signals = cellfun(@(modulation){modulation(:)}, C);
    possible_signals = [possible_signals{:}];

end
