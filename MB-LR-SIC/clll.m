function [H_LR, T] = clll(H, delta)
% CLLL Implements the complex LLL (CLLL) reduction algorithm [1]
%
%   Input  
%       H: Lattice Basis H = [h_1, h_2, ..., h_n] with h_i being a m-dim.
%       column vector
%       delta: Factor to achieve a good quality-complex tradeoff, 1/2 <
%           delta <= 1
%
%   Output
%       H_LR: CLLL-reduced basis
%       T: Unimodular matrix that satifies H_LR = H*T
%
%   [1]: Y. Gan, C. Ling, W. Mow, "Complex Lattice Reduction Algorithm for
%       Low Complexity Full-Diversity MIMO Detection", IEEE Transactions On
%       Signal Processing, Vol. 57, No. 7, July 2009, pp. 2701 - 2710

    if delta <= 1/2
        delta = 0.51;
    elseif delta > 1
            delta = 1;
    end

    n = size(H, 2);  % Iterate over columns in H
    squared_norms = zeros(1, n);  % Denoted with calligraphic H in [1]. The 
                                  % squared norms of the orthoganal vectors   
    for k = 1:n
        squared_norms(k) = dot(H(:,k), H(:,k)); addflops(6 * flops_mul(H(:,k), ... 
            H(:,k)));
    end
    
    % Modified Gram-Schmidt orthogonalization (GSO) procdure to compute the
    % squared norms
    normed_prods = zeros(n,n);  % Denoted with \mu in [1]
    for k = 1:n
        for l = k+1:n 
            inner_prod_kl = dot(H(:,l), H(:,k));
            sub_term = sum(conj(normed_prods(k,1:k-1)) .* ...
                normed_prods(l,1:k-1) .* squared_norms(1:k-1));
            normed_prods(l,k) = (1 / squared_norms(k)) * ... 
                (inner_prod_kl - sub_term);
            squared_norms(l) = squared_norms(l) - ... 
                abs(normed_prods(l, k))^2 * squared_norms(k);
            
            addflops(6 * flops_mul(H(:,l), H(:,k)) + 12 * (k-1) + ... 
                flops_row_sum(1, k-1) + flops_div() + 4 + flops_abs() + ... 
                flops_pow(2) + 2);
        end
    end
    
    T = eye(n);
    k = 2;
    while k <= n
        if abs(real(normed_prods(k,k-1))) > 1/2 || ... 
                abs(imag(normed_prods(k,k-1))) > 1/2
            % First condition (5) for a CLLL-reduced complex lattice is  
            % violated -> Reduce vector size 
            addflops(2 * flops_abs() + 4);
            
            [H, T, normed_prods] = size_reduce(H, T, normed_prods, k, k-1);
        end
        if squared_norms(k) < (delta - abs(normed_prods(k,k-1))^2) * ... 
                squared_norms(k-1)
            % Second condition (6) for a CLLL-reduced complex lattice is
            % violated -> Swap vectors and update ((7) - (15) in [1])
            addflops(3 + flops_abs() + flops_pow(2) + 1);
            
            H(:,[k,k-1]) = H(:,[k-1,k]);  % (7) and (8)
            H_k = squared_norms(k);
            H_k_1 = squared_norms(k-1);
            squared_norms(k-1) = H_k + abs(normed_prods(k,k-1))^2 * ... 
                H_k_1;  % (9)
            normed_prods(k,k-1) = conj(normed_prods(k,k-1)) * ... 
                (H_k_1 / squared_norms(k-1));  % (10)
            squared_norms(k) = H_k_1 - abs(normed_prods(k,k-1))^2 * ... 
                squared_norms(k-1);  % (11)
            
            addflops(1 + flops_abs() + flops_pow(2) + 1 + flops_div() + ... 
                2 + 1 + flops_abs() + flops_pow(2) + 1);
            
            for l = k+1:n
                tmp = normed_prods(l,k-1);
                normed_prods(l,k-1) = normed_prods(l,k-1) * ... 
                    normed_prods(k,k-1) + normed_prods(l,k) * ... 
                    (H_k / squared_norms(k-1));  % (12)
                normed_prods(l,k) = tmp - normed_prods(l,k) * ... 
                    normed_prods(k, k-1);  % (13)
            end
            addflops((n - (k+1)) * 18);
            
            for l = 1:k-2
                normed_prods([k-1,k],l) = normed_prods([k,k-1],l);  % (14)
                % and (15)
            end
            T(:,[k,k-1]) = T(:,[k-1,k]);
            k = max(2,k-1);
        else
            for l = k-2:-1:1
                if abs(real(normed_prods(k,l))) > 1/2 || ... 
                        abs(imag(normed_prods(k,l))) > 1/2
                    % First condition (5) for a CLLL-reduced complex
                    % lattice is violated -> Reduce vector size
                    addflops(2 * flops_abs() + 4);
                    
                    [H, T, normed_prods] = ... 
                        size_reduce(H, T, normed_prods, k, l);
                end
            end
            k = k+1;
        end
    end
    H_LR = H;
    
end