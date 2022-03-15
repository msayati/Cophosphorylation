function [ W ] = mapnetwork( W1, W2, varargin )
    W = logical(double(W1) * double(W2));
    for i = 1:length(varargin)
        W = mapnetwork(W, varargin{i});
    end
end

