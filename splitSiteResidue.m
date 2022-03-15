function [ residue, position ] = splitSiteResidue( A )
    [residue, position] = cellfun(@funSplitResidue, A, ...
        'UniformOutput', false);
    position = cell2mat(position);
end

function [residue, position] = funSplitResidue(x)
    residue = upper(x(1));
%     position = x(regexp(x, '\d'));
    position = str2double(x(regexp(x, '\d')));
end
