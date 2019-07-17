function L = fct_graphLaplacian(A)
    A = (A+A')/2;
    D = diag(sum(A));
    L = D - A;
    D = diag(D) .^ (-0.5);
    D(isinf(D)) = 0;
    D = diag(D);
    L = D*L*D;
end

