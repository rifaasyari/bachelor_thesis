function [H, T, mu] = size_reduce(H, T, mu, k, j)
% SIZE_REDUCE Subroutine in CLLL reduction algorithm [1]. "[...] Aims to
% make basis vectors shorter and closer to orthogonal" [1]. 
% 
%   [1]: Y. Gan, C. Ling, W. Mow, "Complex Lattice Reduction Algorithm for
%       Low Complexity Full-Diversity MIMO Detection", IEEE Transactions On
%       Signal Processing, Vol. 57, No. 7, July 2009, pp. 2701 - 2710
%
%   See also CLLL
    
    c = round(mu(k, j));
    H(:,k) = H(:,k) - c * H(:,j);
    T(:,k) = T(:,k) - c * T(:,j);
    
    addflops(2 * (rows(H) + rows(T)));
    
    for l = 1:j
        mu(k,l) = mu(k,l) - c * mu(j,l);
    end
    
    addflops(l*2);
    
end
