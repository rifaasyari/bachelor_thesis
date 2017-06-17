function [W] = mmse_filter(H)
% MMSE_FILTER Implements the minimum mean-square error (MMSE) filter for
% MIMO systems [1]. This filter minimizes the  error between the transmitted 
% symbol and its estimate at the filter output
%
%   Input
%       H: Complex channel gain matrix
%
%   Output
%       W: MMSE filter which minimizes error between transmitted symbol s
%           and its estimate s_hat = ctranspose(W) * y where y is the 
%           received signal.
%
%   [1]: J. Benesty, Y. Huang, J. Chen, "A Fast Recursive Algorithm for
%   Optimum Sequential Signal Detection in a BLAST System", in IEEE
%   Transactions on Signal Processing, Vol. 51, No. 7, pp. 1722
%   - 1730, July 2003

    W = H*inv(ctranspose(H)*H);  % Calculate the pseudo-inverse of H
    
    addflops(flops_div + flops_mul(rows(H), cols(H), cols(H)) + ...
        flops_inv(cols(H)) + flops_mul(cols(H), rows(H), cols(H)) + ...
        2 * cols(H)^2);
    
end
