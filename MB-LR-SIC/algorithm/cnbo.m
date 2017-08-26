function [P] = cnbo(H, reverse)
% CNBO Implements column-norm based ordering. 
%
%   Input
%       H: Matrix to be ordered
%       reverse (bool): If true sort columns in increasing order of norm.
%           Defaults to false
%
%   Output
%       P: Permutation matrix. H * P results in column-norm based odered 
%           matrix

    if nargin == 1
        reverse = false;
    end
    
    n = size(H,2);  % number of columns in H
    column_norms = zeros(1, n);
    
    for k=1:length(column_norms)
        column_norms(k) = norm(H(:,k));  % Compute norms
        
        addflops(flops_mul(H(:,k)', H(:,k)) + flops_sqrt);
    end
    P = zeros(n, n);
    for k=1:length(column_norms)
        if ~reverse
            [~, arg] = max(column_norms);
            column_norms(arg) = -Inf;  % Ignore index arg for remaining 
                                       % iterations 
        else
            [~, arg] = min(column_norms);
            column_norms(arg) = +Inf;  % Ignore index arg for remaining 
                                       % iterations 
        end
        P(arg, k) = 1;  % Permute column k with column arg

    end  
end
