function [U,V,A,B,cc] = MyCCAPPI(X1, X2, co, perc, L1, L2)
    
    N = size(X2,1);
    
    D1 = size(X1,2);
    D2 = size(X2,2);
    
    if isempty(L1)
        L1 = eye(D1);
    end
    if isempty(L2)
        L2 = eye(D2);
    end
    
    X1 = X1 - repmat(mean(X1),N,1);
    X2 = X2 - repmat(mean(X2),N,1);

    C12 = X1'*X2/N;
    
    C11 = (X1'*X1/N);
    mi = eigs(C11,1,'sr'); if isnan(mi) mi = 0; end
    ma = eigs(C11,1,'lr'); if isnan(mi) ma = 0; end
    lambda = (ma-mi*co)/(co-3);
    C11 = C11 + lambda*L1 + lambda*eye(D1);
    
    C22 = (X2'*X2/N);
    mi = eigs(C22,1,'sr'); if isnan(mi) mi = 0; end
    ma = eigs(C22,1,'lr'); if isnan(mi) ma = 0; end
    lambda = (ma-mi*co)/(co-3);
    C22 = C22 + lambda*L2 + lambda*eye(D2);
    
    A = [zeros(D1,D1) C12; C12' zeros(D2,D2)];
    B = blkdiag(C11, C22);
    clear C11 C22 C12

    [V,D] = eigs(A,B,100,'lr');
    
    d = diag(D);
    d = d(d > 0);
    r = find((cumsum(d) / sum(d)) > perc);
    r = r(1);
    V = V(:,1:r);

    A = V(1:D1,:);
    B = V(D1+1:D1+D2,:);

    
    U = X1*A;
    V = X2*B;
    
    cc = diag(D);
    
end
