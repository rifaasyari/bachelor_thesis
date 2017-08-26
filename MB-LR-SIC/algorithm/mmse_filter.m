function [W] = mmse_filter(H, N, R)
% MMSE_FILTER Implements the minimum mean-square error (MMSE) filter for
% MIMO systems [1]. This filter minimizes the  error between the transmitted 
% symbol and its estimate at the filter output
%
%   Input
%       H: Complex channel gain matrix
%       S: Signal power
%       N: Noise power
%       T_inv: Inverse of uni-modular matrix T used for lattice reduction
%       modulation (str): QPSK or 16-QAM
%       P_l (matrix): Permuation matrix for the l-th branch
%       n (scalar): Index of detected signal component
%       R (matrix): Covariance matrix of symbols in lattice reduced
%           constellation
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
%   [2]: K. Lee, J. Chun, L. Hanzo, "Optimal Lattice-Reduction Aided
%   Successive Interference Cancellation for MIMO Systems", IEEE
%   Transactions on Signal Processing, vol. 6, no. 7, pp. 2438 - 2443, 2007
    
    % P_l = P_l(:, n:end);
    % R = S*(P_l'*(T_inv*T_inv')*P_l);  % R given in [2]

    % W = H * inv(ctranspose(H)*H + ctranspose(N*inv(R(n:end, n:end))));
    W = inv(ctranspose(H)*H+N*inv(R)) * ctranspose(H);
        
%     addflops(flops_div + flops_mul(rows(H), cols(H), cols(H)) + ...
%         flops_inv(cols(H)) + flops_mul(cols(H), rows(H), cols(H)) + ...
%         2 * cols(H)^2);
    
end
