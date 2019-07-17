function r = MyRSquare(y, yhat, m)
    total = sum((y - m) .^ 2);
    residual = sum((y - yhat) .^ 2);
    r = 1 - residual/total;
end

