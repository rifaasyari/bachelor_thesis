function [s_ml] = fsd(y,H,i,constellation,n_S,s_ml,column,l_1)
% FSD Implements the fixed-complexity sphere decoder (FSD). It performs a
%   tree search following down fixed paths and returning the one which
%   satifies the ML criterion.
%
% INPUT y: Received symbol vector. size(y) = N_rx1
%       H: Channel matrix. size(H) = N_rxN_t
%       i: Tree level. Must be initially set to N_t.
%       constellation: Symbol constellation where s, y=Hs+n, is taken from
%       n_S (optional): Distribution of nodes
%       All remaining input arguments are only used for recursive calls of
%           this function
%       
% OUTPUT s_ml: Estimation of transmitted symbol vector s with quasi-ML
%           performance

    P = numel(constellation);
    N_t = size(H,2);
    N_r = length(y);

    if i == N_t
        H_orig = H;
        if nargin == 4
            l_p = ceil(sqrt(N_t)-1);  % Valid iff N_r == N_t
            l_1 = N_t - l_p;
            n_S = [ones(l_1,1); P*ones(l_p,1)];
        else
            l_1 = sum(n_S == 1);
        end
        N_S = prod(n_S);  % Number of searched paths
        s_ml = zeros(N_r,N_S);
        [H,order] = fsd_ordering(H, P, n_S);
        column = 1;
    end

    [~,U] = qr(H);
    s_hat = inv(ctranspose(H)*H)*ctranspose(H)*y;  % Babai point

    s_SE = se_enum(s_hat,s_ml(:,column),U,i,constellation);

    col_end = prod(n_S(i-1:-1:1)) - 1;
    if n_S(i) == P
        for s_i = s_SE
            s_ml(i,column:column+col_end) = s_i;
            if i > 1
                s_ml = fsd(y,H,i-1,constellation,n_S,s_ml,column,l_1);
            end
            column = column + col_end + 1;
        end
    else
        s_ml(i,column:column+col_end) = s_SE(1);
        if i > 1
            s_ml = fsd(y,H,i-1,constellation,n_S,s_ml,column,l_1);
        end
    end

    if i < N_t
        return
    else
        H = H_orig;
        rev_order = arrayfun(@(x)find(order==x,1),1:N_t);
        s_ml = s_ml(rev_order,:);
        Y = repmat(y, 1,N_S);
        ml_arg = Y - H*s_ml;
        ml_criterion = dot(ml_arg,ml_arg);
        [~, arg_min] = min(ml_criterion);
        s_ml = s_ml(:,arg_min);
    end

end