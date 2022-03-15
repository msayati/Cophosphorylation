function [W] = generatePositionMatrix( n, dist )
    W = false(n, n);
    for i = 1:n
        iMax = min(i + dist, n);
        W(i, i:iMax) = true;
    end
end

