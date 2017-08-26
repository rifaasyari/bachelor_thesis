function [Q,R] = psa(Q,R,N_t,N_r)
    
    Q = Q(:,1:4);
    R = R(1:4,:);

    Q1 = Q(1:N_r,:);
    Q2 = Q(N_r+1:end,:);
    
    for i = N_t:-1:2
        row_norms = sqrt(dot(Q2(1:i,1:i),Q2(1:i,1:i)));
        Q2([i,argmin(row_norms)],:) = Q2([argmin(row_norms),i],:);
        R([i,argmin(row_norms)],:) = R([argmin(row_norms),i],:);
        
        v = R(i,:);
        Theta = eye(N_t) - 2*(v*v');
        
        R = R*Theta;
        Q1 = Q1*Theta;
    end
    Q = [Q1,Q2];
end