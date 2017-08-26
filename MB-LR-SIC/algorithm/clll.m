function [H_LR, T, P_c] = clll(H, delta)
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
%       P_c: Probability that test in line 107 is passed 
%
%   [1]: Y. Gan, C. Ling, W. Mow, "Complex Lattice Reduction Algorithm for
%       Low Complexity Full-Diversity MIMO Detection", IEEE Transactions On
%       Signal Processing, Vol. 57, No. 7, July 2009, pp. 2701 - 2710

    test_arrival = 0;
    test_passed = 0;

    if delta <= 1/2
        delta = 0.51;
    elseif delta > 1
            delta = 1;
    end
    
    n = size(H, 2);  % Iterate over columns in H
    squared_norms = zeros(1, n);  % Denoted with calligraphic H in [1]. The 
                                  % squared norms of the original vectors   
    for j = 1:n
        squared_norms(j) = dot(H(:,j), H(:,j)); 
        
        addflops(flops_mul(H(:,j)',H(:,j)));
        
    end
    
    % Modified Gram-Schmidt orthogonalization (GSO) procdure to compute the
    % squared norms
    mu = ones(n,n);  % Denoted with \mu in [1]. Initialize with zeros or 
                     % with ones?
    for j = 1:n
        for i = j+1:n 
            inner_prod_ij = dot(H(:,j), H(:,i));
            sub_term = sum(conj(mu(j,1:j-1)) .* ...
                mu(i,1:j-1) .* squared_norms(1:j-1));
            mu(i,j) = (1 / squared_norms(j)) * ... 
                (inner_prod_ij - sub_term);
            squared_norms(i) = squared_norms(i) - ... 
                abs(mu(i, j))^2 * squared_norms(j);
            
            addflops(flops_mul(H(:,i)', H(:,j)) + 2 * (j-1) + ... 
                flops_row_sum(1, j-1) + flops_div() + 4 + flops_abs() + ... 
                flops_pow(2) + 2);
        end
    end
    
    T = eye(n);
    k = 2;
    while k <= n
        if abs(real(mu(k,k-1))) > 1/2 || abs(imag(mu(k,k-1))) > 1/2
            % First condition (5) for a CLLL-reduced complex lattice is  
            % violated -> Reduce vector size 
            addflops(2 * flops_abs() + 4);
            
            [H, T, mu] = size_reduce(H, T, mu, k, k-1);
        end
        if squared_norms(k) < (delta - abs(mu(k,k-1))^2) * ... 
                squared_norms(k-1)
            % Second condition (6) for a CLLL-reduced complex lattice is
            % violated -> Swap vectors and update ((7) - (15) in [1])
            addflops(3 + flops_abs() + flops_pow(2) + 1);
            
            H(:,[k,k-1]) = H(:,[k-1,k]);  % (7) and (8)
            H_k = squared_norms(k);
            H_k_1 = squared_norms(k-1);
            squared_norms(k-1) = H_k + abs(mu(k,k-1))^2 * ... 
                H_k_1;  % (9)
            mu_old = mu;
            mu(k,k-1) = conj(mu(k,k-1)) * ... 
                (H_k_1 / squared_norms(k-1));  % (10)
            squared_norms(k) = H_k_1 - abs(mu(k,k-1))^2 * ... 
                squared_norms(k-1);  % (11)
            
            addflops(1 + flops_abs() + flops_pow(2) + 1 + flops_div() + ... 
                2 + 1 + flops_abs() + flops_pow(2) + 1);
            
            for i = k+1:n
                mu(i,k-1) = mu_old(i,k-1) * ... 
                    mu(k,k-1) + mu_old(i,k) * ... 
                    (H_k / squared_norms(k-1));  % (12)
                mu(i,k) = mu_old(i,k-1) - mu_old(i,k) * ... 
                    mu_old(k, k-1);  % (13)
                
                addflops(flops_div() + 5);
            end
            
            for j = 1:k-2
                mu(k-1,j) = mu_old(k,j);  % (14) 
                mu(k,j) = mu_old(k-1,j);  % (15)
            end
            T(:,[k,k-1]) = T(:,[k-1,k]);
            k = max(2,k-1);
        else
            for j = k-2:-1:1
                test_arrival = test_arrival + 1;
                if abs(real(mu(k,j))) > 1/2 || abs(imag(mu(k,j))) > 1/2
                    % First condition (5) for a CLLL-reduced complex
                    % lattice is violated -> Reduce vector size
                    addflops(2 * flops_abs() + 4);
                    test_passed = test_passed + 1;
                    
                    [H, T, mu] = size_reduce(H, T, mu, k, j);
                end
            end
            k = k+1;
        end
    end
    H_LR = H;
    P_c = test_passed / test_arrival;  % Probability that conditional test 
                                       % is passed
    
end