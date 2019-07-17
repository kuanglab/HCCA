function [U,V,W,A,B,C] = MyCCA3(X1, X2, X3, co, perc)
    D1 = size(X1,2);
    D2 = size(X2,2);
    D3 = size(X3,2);
    N = size(X2,1);
    
    X1 = X1 - repmat(mean(X1),N,1);
    X2 = X2 - repmat(mean(X2),N,1);
    X3 = X3 - repmat(mean(X3),N,1);
    
    C12 = X1'*X2/N;
    C13 = X1'*X3/N;
    C23 = X2'*X3/N;
    
    C11 = (X1'*X1/N);
    mi = eigs(C11,1,'sr'); if isnan(mi) mi = 0; end
    ma = eigs(C11,1,'lr'); if isnan(mi) ma = 0; end
    lambda = (ma-mi*co)/(co-1);
    C11 = C11 + lambda*eye(size(X1,2));
    
    C22 = (X2'*X2/N);
    mi = eigs(C22,1,'sr'); if isnan(mi) mi = 0; end
    ma = eigs(C22,1,'lr'); if isnan(mi) ma = 0; end
    lambda = (ma-mi*co)/(co-1);
    C22 = C22 + lambda*eye(size(X2,2));
    
    C33 = (X3'*X3/N);
    mi = eigs(C33,1,'sr'); if isnan(mi) mi = 0; end
    ma = eigs(C33,1,'lr'); if isnan(mi) ma = 0; end
    lambda = (ma-mi*co)/(co-1);
    C33 = C33 + lambda*eye(size(X3,2));

    A = [C11 C12 C13; C12' C22 C23; C13' C23' C33];
    B = blkdiag(C11, C22, C33);
%     [V,D] = eig(A,B);
%     [~,ind] = sort(real(diag(D)));
%     V = V(:,ind(end:-1:end-r+1));

    [V,D] = eigs(A,B,100,'lr');
    D = D - eye(100);
    
    d = diag(D);
    r = find((cumsum(diag(D)) / sum(diag(D))) > perc);
    r = r(1);
    V = V(:,1:r);
    
    A = V(1:D1,:);
    B = V(D1+1:D1+D2,:);
    C = V(D1+D2+1:D1+D2+D3,:);
    
    U = X1*A;
    V = X2*B;
    W = X3*C;
    
end

