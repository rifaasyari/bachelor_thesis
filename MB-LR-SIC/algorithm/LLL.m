function [Q_tilde, R_tilde, T] = LLL(Q,R, delta, P)

    m = size(R,2);  % m = 2*N_t

    if nargin == 3
        P = eye(m);
    end
    
    Q_tilde = Q;
    R_tilde = R;
    T = P;
    
    k = 2;
    
    while k <= m
        
        for l = k-1:-1:1
            mu = round(R_tilde(l,k)/R_tilde(l,l));
            if mu ~= 0
                R_tilde(1:l,k) = R_tilde(1:l,k) - mu*R_tilde(1:l,l);
                T(:,k) = T(:,k) - mu*T(:,l);
            end
        end
        
        if delta*R_tilde(k-1,k-1)^2 > R_tilde(k,k)^2 + R_tilde(k-1,k)^2
            R_tilde(:,[k-1,k]) = R_tilde(:,[k,k-1]);
            T(:,[k-1,k]) = T(:,[k,k-1]);
            alpha = R_tilde(k-1,k-1)/norm(R_tilde(k-1:k,k-1));
            beta = R_tilde(k,k-1)/norm(R_tilde(k-1:k,k-1));
            Theta = [alpha beta; -beta alpha];
            R_tilde(k-1:k,k-1:m) = Theta*R_tilde(k-1:k,k-1:m);
            Q_tilde(:,k-1:k) = Q_tilde(:,k-1:k)*Theta';
            k = max(k-1,2);
        else
            k = k+1;
        end
        
    end

end