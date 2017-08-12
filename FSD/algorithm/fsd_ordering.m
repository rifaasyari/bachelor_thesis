function [H_o, order] = fsd_ordering(H,P,n_S)
% FSD_ORDERING Implements the FSD ordering of the channel matrix. Ordering
%   is done in a way that the signal with the highest post-detection
%   noise amplification is selected in level i if n_S(i) == P (where P
%   refers to the number of symbols in the used symbol constellation for
%   transmission). Otherwise the signal with the lowest post-deteciton
%   noise amplification is selected.

% INPUT H: Channel matrix. size(H) = N_rxN_t
%       P: Number of constellation points
%       n_S: Distribution of nodes
%
% OUTPUT H_o: FSD ordered channel matrix

    % TODO: Look for more details on FSD ordering

    H_o = zeros(size(H));
    N_t = size(H,2);
    order = size(1,N_t);
    cols = 1:N_t;
    
    for i = N_t:-1:1
        H_inv = inv(ctranspose(H)*H)*ctranspose(H);
        noise_amplif = dot(H_inv,H_inv,2);
        if n_S(i) == P
            [~, k] = max(noise_amplif);
        else
            [~, k] = min(noise_amplif);
        end
        H_o(:,i) = H(:,k);
        order(i) = cols(k);
        H = [H(:,1:k-1), H(:,k+1:end)];
        cols(k) = [];
    end

end