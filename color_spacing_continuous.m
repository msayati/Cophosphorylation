function [C] = color_spacing_continuous( values, breaks, colors)
    indices = zeros(size(values));
    for i = 1:length(breaks)
       if(i ~= length(breaks))
           q = values >= breaks(i);
       else
           q = values > breaks(i);
       end
       indices(q) = indices(q) + 1;
    end
    max(indices)
    ratio = zeros(size(values));
    ranges = breaks(2:end) - breaks(1:(end)-1);
    for i = 1:(length(breaks)-1)
        ind = indices == i;
        ratio(ind) = (values(ind) - breaks(i)) / (ranges(i));        
    end
    C = (1-ratio).*colors(indices, :) + ratio.*colors(indices+1, :);
end

