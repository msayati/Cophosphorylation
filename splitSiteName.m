function [ protein, site ] = splitSiteName( A, delimiter, joinFirst)
    if(nargin < 3)
        joinFirst = true;
    end
    if(joinFirst)
        [protein, site] = cellfun(@(x) funSplit(x, delimiter), A, ...
            'UniformOutput', false);
    else
        [protein, site] = cellfun(@(x) funSplit2(x, delimiter), A, ...
            'UniformOutput', false);
    end
end

function [protein, site] = funSplit(x, delimiter)
    c = strsplit(x, delimiter);
    if(length(c) == 1)
       error(sprintf('Delimited not found: %s', x));
    end
    if(length(c) > 2)
        protein = char(join(c(1:(end-1)), delimiter));
        site = c{end};
    else
        protein = c{1};
        site = c{2};
    end
end

function [protein, site] = funSplit2(x, delimiter)
    c = strsplit(x, delimiter);
    if(length(c) == 1)
       error(sprintf('Delimited not found: %s', x));
    end
    if(length(c) > 2)
        protein = c{1};
        site = char(join(c(2:(end)), delimiter));
    else
        protein = c{1};
        site = c{2};
    end
end
