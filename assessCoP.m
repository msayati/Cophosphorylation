function [D] = assessCoP(X, W, measure)
    A = cophospho(X, 'measure', measure);
    m = size(W, 1);
    [i, j] = ind2sub(size(W), find(triu(W, 1)));
    D = A((i-1).*(m-i/2)+j-i);
end

