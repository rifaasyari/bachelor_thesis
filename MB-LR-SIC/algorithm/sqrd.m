function [P,Q,R] = sqrd(H)

    R_square = H'*H;
    diag_R = diag(R_square);
    [~,column_order] = sort(diag_R);

    P = eye(size(H,2));
    P = P(:,column_order);
    
    [Q,R] = qr(H*P);
end
