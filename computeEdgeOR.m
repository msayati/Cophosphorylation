function [LOR, LOR_err] = computeEdgeOR( labels )
    nEdge = length(labels);
    nLabel = nnz(labels);
    
    cum_labels = cumsum(labels);
    num_edges = cumsum(ones(size(labels)));
    TT = cum_labels;
    TF = num_edges - cum_labels;
    FT = nLabel - cum_labels;
    FF = nEdge - nLabel - num_edges + cum_labels;
    OR = (FF .* TT) ./ (TF .* FT);
    LOR = log2(OR);
    LOR_err = 1.96 * sqrt(1./FF + 1./TT + 1./TF + 1./FT) / log(2);
    invalids = (TT <= 0) | (TF <= 0) | (FT <= 0) | (FF <= 0);
    LOR(invalids) = NaN;
    LOR_err(invalids) = NaN;
end

