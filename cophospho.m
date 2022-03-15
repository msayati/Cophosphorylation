function [ D ] = cophospho( X, varargin )
    p = inputParser;
    p.CaseSensitive = false;
    validMatrix = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'2d', 'nonempty', 'real', 'nonsparse', 'nonnan', 'finite'});
    validMeasure = @(x) any(validatestring(x, ...
        {'absdif', 'abscorr', 'corr', 'bicorr'}));
    addRequired(p, 'X', validMatrix);
    addParameter(p, 'Measure', 'absdif', validMeasure);
    parse(p, X, varargin{:});
    param = p.Results;
    
    nSite = size(X, 1);
    nSample = size(X, 2);
    
    switch(lower(param.Measure))
        case 'corr'
            D = 1 - pdist(X, 'correlation');
        case 'abscorr'
            D = abs(1 - pdist(X, 'correlation'));
        case 'bicorr'
            D = Biweight_midcorrelation(X);
            D = D - diag(diag(D));
            D = squareform(D);
        case 'absdif'
            D = pdist(X, 'cityblock') ./ nSample;
        otherwise
            error('Invalid co-phosphorylation measure.');
    end
    
end

