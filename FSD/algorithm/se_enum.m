function nodes = se_enum(s_hat,s,U,i,nodes)
% SE_ENUM Performs the Schnorr-Euchner (SE) enumeration for nodes on level
%   i. The SE enumeration orders nodes in increasing distance to the 
%   hypersphere radius z_i on level i.
%
% INPUT s_hat: Babai point
%       s: Nodes on the current path
%       U: N_txN_t-dim. upper triangular matrix for QR decomposition of H
%       i: Current tree level
%       nodes: All children on level i
%
% OUTPUT nodes: Child nodes in SE enumeration

    % Setup    
    S = repmat(s, 1, length(nodes));
    S(i,:) = nodes';
    S_hat = repmat(s_hat, 1, length(nodes));
    
    % Accumulated squared euclidean distance
    P_i = abs(U(i,i:end)*(S(i:end,:)-S_hat(i:end,:))).^2;
    
%     addflops(2*numel(S(i:end,:)) + flops_mul(U(i,i:end),S(i:end,:)) + ...
%         flops_abs() * numel(nodes) + numel(nodes)*flops_pow(2));
    
    % Sort according to increasing euclidean distance 
    [~, order] = sort(P_i, 'ascend');
    nodes = nodes(order)';
end
