function [S_ml, R, E_n, visit_children, full_paths, visited] = ... 
    sd(y,H,i,constellation,R,P_prev,S_ml,E_n, visit_children, s)
% SD Implementing the complex version of the Sphere Decoder for MIMO
%   detection. Achieves quasi-ML performance under variable complexity.
%
% INPUT y: Received symbol vector
%       H: N_rxN_t-dim. channel gain matrix
%       i: Current tree level. Initially: i=N_t
%       constellation: Symbol constellation used for transmission
%       The remaining input arguments are used for recursive calls of this
%           function
%
% OUTPUT s_ml: ML solution to y
%        R: Updated hypersphere radius corresponding to s_best
%        E_n: Average number of child nodes visited per parent node per 
%           level i from i=N_t to i=1
%        full_paths: Number of paths visited to the leaf nodes
%        visited: Number of paths which were not fully visited

    % Tree setup
    N_t = size(H,2);

    [~, U] = qr(H);
    s_hat = inv(ctranspose(H)*H)*ctranspose(H)*y;  % Babai point
        
    if i == N_t
        E_n = zeros(N_t,1);  % Average number of nodes per level 
        visit_children = [zeros(N_t-1,1);1];  % How many parent nodes visit 
        % at least one child node?
        
        R=Inf;
        P_prev=0;
        S_ml = {};  % Store all possible ML paths
        s = zeros(N_t,1);
    end
    
    s_SE = se_enum(s_hat, s, U, i, constellation);  
    % SE enumeration of nodes
    
    first_child = true;
    
    % Depth-first search
    for s_i = s_SE
        
        if i == N_t
            P_i = abs((s_i-s_hat(N_t)))^2*U(N_t,N_t)^2;  
            % Path metric on level i=N_t
        else
            P_i = abs(U(i,i:end)*([s_i;s(i+1:end)]-s_hat(i:end)))^2 ... 
                + P_prev;  
            % Path metric on level i < N_t 
        end
        
        if P_i <= R^2  % Path lies inside the hypersphere
            s(i) = s_i;
            E_n(i) = E_n(i) + 1;
            if i < N_t && first_child
                visit_children(i) = visit_children(i) + 1;
                first_child = false;
            end
            if i == 1  
                % Last level reached, return with ML path and updated
                % hypersphere radius
                % R = norm(U*(S_ml(:,col)-s_hat));
                S_ml{end+1} = s;
                R = sqrt(P_i);
            else
                [S_ml, R, E_n, visit_children] = ... 
                    sd(y,H,i-1,constellation,R,P_i,S_ml,E_n,visit_children,s);
                % Recursively visiting all child nodes, thus realizing 
                % depth-first search
            end
        end
    end
    
    if i < N_t
        % All children visited. Go back to i+1, i.e. one level up
        return
    else
        % Return ML solution
        S_ml = cell2mat(S_ml);
        visited = size(S_ml(:, any(S_ml)),2);
        S = S_ml(:, all(S_ml));  % Discard columns with zero symbols
        S_ml = S(:,end);  % Last path found inside the hypersphere 
                          % corresponds to the ML solution
                          
        full_paths = size(S,2);
        E_n = flip(E_n ./ visit_children);
    end
end
