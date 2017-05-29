function [P] = cnbo(H)
% CNBO Implements column-norm based ordering. 
%
%   Input
%       H: Matrix to be ordered
%
%   Output
%       P: Permutation matrix. H * P results in column-norm odered matrix

    n = size(H,2);  % number of columns in H
    column_norms = zeros(1, n);
    for k=1:length(column_norms)
        column_norms(k) = norm(H(:,k));  % Compute norms
    end
    P = zeros(n, n); 
    for k=1:length(column_norms)
        [argvalue, argmax] = max(column_norms);
        P(argmax, k) = 1;  % Permutate column k with column argmax
        column_norms(argmax) = 0;  % Ignore max_index for remaining 
                                   % iterations 
    end  
end
