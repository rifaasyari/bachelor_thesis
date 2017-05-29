function [W] = mmse_filter(H, signal_power, noise_power)
% MMSE_FILTER Implements the minimum mean-square error (MMSE) filter for
% MIMO systems [1]. This filter minimizes the  error between the transmitted 
% symbol and its estimate at the filter output
%
%   Input
%       H: Complex channel gain matrix
%       signal_power: Power of transmitted signals
%       noise_power: Power of noise over the channel
%
%   Output
%       W: MMSE filter which minimizes error between transmitted symbol s
%           and its estimate ŝ = ctranspose(W) * y where y is the received
%           signal.
%
%   [1]: J. Benesty, Y. Huang, J. Chen, "A Fast Recursive Algorithm for
%   Optimum Sequential Signal Detection in a BLAST System", in IEEE
%   Transactions on Signal Processing, Vol. 51, No. 7, pp. 1722
%   - 1730, July 2003

    alpha = noise_power / signal_power;  % noise-to-signal-power ratio
    W = H * inv((ctranspose(H) * H + alpha * eye(size(H, 2))));
    
    addflops(flops_div + 2 * flops_mul(rows(H), cols(H), cols(H)) + ...
        flops_inv(cols(H)) + 6 * flops_mul(cols(H), rows(H), cols(H)) + ...
        2 * cols(H)^2);
end
