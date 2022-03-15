function [mValues] = assessPairwiseMean(V, W, measure)
    if(nargin < 3)
       measure = 'absdif'; 
    end
    switch(measure)
        case 'absdif'
            measure_id = 1;
        case 'corr'
            measure_id = 2;
        otherwise
            error('Invalid co-phosphorylation measure.');
    end

    [p1, p2] = ind2sub(size(W), find(W));
    if(issymmetric(W))
        valids = p1 < p2;
        p1 = p1(valids);
        p2 = p2(valids);
    end
    mValues = zeros(length(p1), 1);
    for iPair = 1:length(p1)
        i1 = p1(iPair);
        i2 = p2(iPair);
        p1v = V(i1, :);
        p2v = V(i2, :);
        if(measure_id == 1)
            mValues(iPair) = mean(abs(p1v - p2v), 'omitnan');
        end
        if(measure_id == 2)
            mValues(iPair) = corr(p1v', p2v', 'rows','complete');
        end
    end
%     mValues = mValues(~isnan(mValues));
end