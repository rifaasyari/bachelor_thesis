function P = pno(H)
% PNO Post Noise detection Ordering.
%   Ordering is done in a way that the signal with the lowest post-detection
%   noise amplification are detected first

% INPUT H: Channel matrix. size(H) = N_rxN_t
%       P: Number of constellation points
%       n_S: Distribution of nodes
%
% OUTPUT H_o: FSD ordered channel matrix

    H_o = zeros(size(H));
    N_t = size(H,2);
    P = zeros(N_t);
    cols = 1:N_t;
    
    for i = 1:N_t
        H_inv = inv(ctranspose(H)*H)*ctranspose(H);
        noise_amplif = dot(H_inv,H_inv,2);
        [~, k] = min(noise_amplif);
        H_o(:,i) = H(:,k);
        P(cols(k),i) = 1;
        H = [H(:,1:k-1), H(:,k+1:end)];
        cols(k) = [];
    end

end