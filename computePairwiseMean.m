function [mValues] = computePairwiseMean(W, p1, p2)
	if(nargin < 2)
		[p1, p2] = ind2sub(size(W), find(W));
		if(issymmetric(W))
			valids = p1 < p2;
			p1 = p1(valids);
			p2 = p2(valids);
		end
	end
    mValues = zeros(length(p1), 1);
    for iPair = 1:length(p1)
%         if(mod(iPair, 100000) == 0)
%            disp(iPair); 
%         end
        i1 = p1(iPair);
        i2 = p2(iPair);
        p1v = W(i1, :);
        p2v = W(i2, :);
        mValues(iPair) = mean(abs(p1v - p2v), 'omitnan');
    end
    mValues = mValues(~isnan(mValues));
end