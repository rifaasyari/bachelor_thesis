function [H, T, normed_prods] = size_reduce(H, T, normed_prods, k, l)
% SIZE_REDUCE Subroutine in CLLL reduction algorithm [1]. "[...] Aims to
% make basis vectors shorter and closer to orthogonal" [1]. 
% 
%   [1]: Y. Gan, C. Ling, W. Mow, "Complex Lattice Reduction Algorithm for
%       Low Complexity Full-Diversity MIMO Detection", IEEE Transactions On
%       Signal Processing, Vol. 57, No. 7, July 2009, pp. 2701 - 2710
%
%   See also CLLL

    c = round(normed_prods(k, l));
    H(:,k) = H(:,k) - c * H(:,l);
    T(:,k) = T(:,k) - c * U(:,l);
    for m = 1:l
        normed_prods(k,m) = normed_prods(k,m) - c * normed_prods(l,m);
    end
