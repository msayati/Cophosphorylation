function [hP, hS] = violinplot( X, varargin )


    p = inputParser;
    validX = @(x) validateattributes(x, {'cell', 'numeric', 'logical'}, ...
        {'2d', 'nonempty', 'real'});
    validNormalizedScalar = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar','nonempty','real','nonnan','positive', '<=', 1});
    validPosScalar = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar','nonempty','real','nonnan','positive'});
    validColorArray = @(x) validateattributes(x, {'numeric'}, ...
        {'2d', 'nonempty', 'real', 'nonnan', 'ncols', 3});
%     validColor = @(x) validateattributes(x, {'numeric'}, {'vector', ...
%         'nonempty', 'real', 'nonnan', 'ncols', 3, 'nrows', 1, '<=', 1});
    validNamingArray = @(x) validateattributes(x, {'cell'}, ...
        {'vector', 'nonempty'});
%     validStyle = @(x) validateattributes(x, {'char'}, ...
%         {'vector', 'nonempty'});
    validFront = @(x) any(validatestring(x, {'violin', 'scatter'}));
    validPosition = @(x) ...
        any(validatestring(x, {'left', 'center', 'right'}));
    validLineStyle = @(x) any(validatestring(x, {'-', '--', '-.', ':'}));
    
    addRequired(p, 'X', validX);
    addParameter(p, 'XShift', 0, @isnumeric);
    addParameter(p, 'Clipped', true, @islogical);
    addParameter(p, 'Darkening', 0.5, validNormalizedScalar);
    addParameter(p, 'FontSize', 12, validPosScalar);
    addParameter(p, 'Jitter', 0.5, @validateJitter);
    addParameter(p, 'JitterPosition', 'center', validPosition);
    addParameter(p, 'Front', 'scatter', validFront);
    addParameter(p, 'Group', [], validX);
    addParameter(p, 'GroupColors', [0.3 0.6 0.8], validColorArray);
    addParameter(p, 'GroupNames', [], validNamingArray);
    addParameter(p, 'LineWidth', 0.75, validPosScalar);
    addParameter(p, 'MarkerAlpha', 1, validNormalizedScalar);
    addParameter(p, 'MarkerColor', [0.1 0.1 0.1], validColorArray);
    addParameter(p, 'MarkerLineWidth', 1, validPosScalar);
    addParameter(p, 'MarkerSize', 7, validPosScalar);
%     addParameter(p, 'MarkerStyle', '+', validStyle);
    addParameter(p, 'ViolinAlpha', 1, validNormalizedScalar);
    addParameter(p, 'ViolinLineColor', 'darken', @validateLineColor);
    addParameter(p, 'ViolinPosition', 'center', validPosition);
    addParameter(p, 'Violin', true, @islogical);
    addParameter(p, 'NormCoef', 'auto', @validateNormCoef);
    addParameter(p, 'Scatter', false, @islogical);
    addParameter(p, 'ShowMean', false, @islogical);
    addParameter(p, 'Boundaries', 'unbounded', @(x) true);
    addParameter(p, 'MeanLine', 'dynamic', @validateJitter);
    addParameter(p, 'MeanLineStyle', '--', validLineStyle);
    addParameter(p, 'MeanLineWidth', 1.25, validPosScalar);
    addParameter(p, 'MeanLineColor', 'darken', @validateLineColor);
    addParameter(p, 'Width', 0.8, validNormalizedScalar);
    parse(p, X, varargin{:});
    param = p.Results;
    checkUsingDefaults = @(p,varname) any(strcmp(p.UsingDefaults,varname));
    
    if(iscell(X))
        nObs = 0;
        for iGroup = 1:length(X)
            xx = X{iGroup};
            nObs = nObs + numel(xx);
            X{iGroup} = reshape(xx, numel(xx), 1);
        end
        
        Xt = nan(nObs, 1);
        param.Group = nan(nObs, 1);
        current = 1;
        for iGroup = 1:length(X)
            xx = X{iGroup};
            indices = current:(current + numel(xx) - 1);
            Xt(indices) = xx;
            param.Group(indices) = iGroup;
            current = current + numel(xx);
        end
        groupsAreSet = true;
        X = Xt;
        clear Xt
    else
        groupsAreSet = false;
    end
    
    if(isempty(X))
       warning('Ignoring empty matrix. Returning without action.');
       hP = [];
       hS = [];
       return;
    end
    
    if(~groupsAreSet && checkUsingDefaults(p, 'Group'))
        [nObs, nGroup] = size(X);
        param.Group = repmat(1:nGroup, nObs, 1);
        param.Group = param.Group(:);
    else
        [GroupNames, ~, param.Group] = unique(param.Group);
        nGroup = length(GroupNames);
    end
    
    if(checkUsingDefaults(p, 'GroupNames'))
        param.GroupNames = cellfun(@(str) num2str(str), ...
            num2cell(1:nGroup), 'UniformOutput', false);
    end
    param.NumGroup = nGroup;
    
    checkErrors(X, param);
    
    if(size(param.GroupColors, 1) == 1)
        param.GroupColors = repmat(param.GroupColors, nGroup, 1);
    end
    
    if(size(param.MarkerColor, 1) == 1)
        param.MarkerColor = repmat(param.MarkerColor, nGroup, 1);
    end
    
    hP = gobjects(nGroup, 1);
    if(param.Scatter)
       hS = gobjects(nGroup, 1);
    else
       hS = gobjects(0, 1);
    end
    
    if(ishold())
        holdFlag = false;
    else
        cla();
        hold('on');
        holdFlag = true;
    end
    
    fMaxes = nan(1, nGroup);
    for iGroup = 1:nGroup
        ind = param.Group == iGroup;
        Xv = X(ind);
        [f] = runKSdensity(Xv, param);
        fMax = max(f);
        fMaxes(iGroup) = fMax;
    end
    
    if(~ischar(param.NormCoef))
       if(length(param.NormCoef) == 1)
           param.NormCoef = repmat(param.NormCoef, nGroup, 1);
           fMaxes = param.NormCoef;
       else
           if(length(param.NormCoef) == nGroup)
               % Check for errors ....
               fMaxes = param.NormCoef;
           else
               error(['Number of normalization coefficients must be ', ...
                   'equal to the number of groups.']);
           end
       end
    end
    
    for iGroup = 1:nGroup
        ind = param.Group == iGroup;
        Xv = X(ind);
%         [f, xi] = ksdensity(Xv);
        [f, xi] = runKSdensity(Xv, param);
        xMid = param.XShift + iGroup - 0.5;
        F = (f / fMaxes(iGroup)) * param.Width / 2;
        Xo = xi;
        Fo = F;
        if(param.Clipped)
            F = [0, F, 0];
            xi = [min(Xv), xi, max(Xv)];
            xi(xi < min(Xv)) = min(Xv);
            xi(xi > max(Xv)) = max(Xv);
        end
        fLeft = xMid - F;
        fRight = xMid + F;
%         Xv = [min(Xv); Xv; max(Xv)];
        
%         xJitter = xMid + (rand(size(Xv)) - 0.5) * param.Jitter;
        
        xJitter = xMid + computeJitter(param, Fo, Xo, Xv);
        faceColor = param.GroupColors(iGroup, :);
        
        patchOpts = struct();
        patchOpts.LineWidth = param.LineWidth;
        patchOpts.EdgeColor = getViolinLineColor(param, faceColor);
        
        plotOpts = struct();
        plotOpts.LineWidth = param.MeanLineWidth;
        plotOpts.Color = getViolinLineColor(param, faceColor);
        
        switch(lower(param.ViolinPosition))
            case 'left'
                Vx = fLeft;
                Vy = xi;
            case 'center'
                Vx = [fLeft flip(fRight)];
                Vy = [xi flip(xi)];
            case 'right'
                Vx = flip(fRight);
                Vy = flip(xi);
            otherwise
                error('Invalid violin position.');
        end
%         c = reshape(faceColor, 1, 1, 3);
%         C = repmat(c, size(Vx, 1), size(Vx, 2));
%         C(1, 1, :) = [1 0 0];
%         C(1, end, :) = [1 0 0];
%         size(Vx)
        if(param.Violin && ~strcmpi(param.Front, 'violin'))
            hPatch = patch(Vx, Vy, faceColor, patchOpts);
        end
        
        if(param.Scatter)
            markerColor = param.MarkerColor(iGroup, :);
            hScatter = scatter(xJitter, Xv, param.MarkerSize, ...
                markerColor, 'filled');
            alpha(hScatter, param.MarkerAlpha);
            hS(iGroup) = hScatter;
        end
        
        if(param.Violin && strcmpi(param.Front, 'violin'))
            hPatch = patch(Vx, Vy, faceColor, patchOpts);
        end
        if(param.Violin)
            alpha(hPatch, param.ViolinAlpha);
            hP(iGroup) = hPatch;
        end
        
        if(param.ShowMean)
            [m, xJ] = computeMean(param, F, xi, Xv);
            plotOpts.Color = getMeanLineColor(param, faceColor);
            plot(xMid + xJ, [m, m], param.MeanLineStyle, plotOpts);
        end
    end
    if(holdFlag)
        hold('off');
    end
    xlim([0 nGroup]);
    set(gca, 'XTick', (1:nGroup) - 0.5, 'XTickLabel', param.GroupNames);
end

function [f, xi] = runKSdensity(Xv, param)
    [f, xi] = ksdensity(Xv, 'Support', param.Boundaries);
%     [f, xi] = ksdensity(Xv, 'Function', 'cdf');
end

function [] = checkErrors(X, param)
    if(numel(param.Group) ~= numel(X))
        error(['Number of elements in grouping variable must be ' ...
            'equal to the number of elements in X.']);
    end
    
    if(numel(param.GroupNames) ~= param.NumGroup)
        error(['Length of group names array must be ' ...
            'equal to the number of groups.']);
    end
    
end

function [color] = getMeanLineColor(param, groupColor)
    color = getLineColor(param.MeanLineColor, groupColor, param.Darkening);
end

function [color] = getViolinLineColor(param, groupColor)
    color = getLineColor(param.ViolinLineColor, groupColor, param.Darkening);
end

function [color] = getLineColor(colorParam, groupColor, darkening)
    if(ischar(colorParam))
        switch(lower(colorParam))
            case 'darken'
                color = brighten(groupColor, -1 * darkening);
            case 'groupcolor'
                color = groupColor;
            otherwise
                error('Invalid line coloring option.');
        end
    else
        color = colorParam;
    end
end

function [] = validateLineColor(x)
    if(ischar(x))
        validatestring(x, {'darken', 'groupcolor'});
    else
        validateattributes(x, {'numeric'}, ...
            {'2d', 'nonempty', 'real', 'nonnan', 'ncols', 3})
    end
end

function [m, xJ] = computeMean(param, F, Xi, Xv)
    switch(lower(param.ViolinPosition))
        case 'left'
            meanLinePos = -0.5 * [0, 1];
        case 'center'
            meanLinePos = ([0, 1] - 0.5);
        case 'right'
            meanLinePos = +0.5 * [0, 1];
        otherwise
            error('Invalid jitter position.'); 
    end
    m = mean(Xv);
    if(ischar(param.MeanLine))
        switch(lower(param.MeanLine))
            case 'dynamic'
                Ix = knnsearch(Xi', m);
                xJ = 2 * meanLinePos .* F(Ix)';
            otherwise
                error('Invalid mean line option.');
        end
    else
        xJ = meanLinePos * param.MeanLine;
    end
end


function [xJ] = computeJitter(param, F, Xi, Xv)
    switch(lower(param.JitterPosition))
        case 'left'
            jitterRNG = -0.5 * rand(size(Xv));
        case 'center'
            jitterRNG = (rand(size(Xv)) - 0.5);
        case 'right'
            jitterRNG = +0.5 * rand(size(Xv));
        otherwise
            error('Invalid jitter position.'); 
    end

    if(ischar(param.Jitter))
        switch(lower(param.Jitter))
            case 'dynamic'
                Ix = knnsearch(Xi', Xv);
                xJ = 2 * jitterRNG .* F(Ix)';
            otherwise
                error('Invalid jitter option.');
        end
    else
        xJ = jitterRNG * param.Jitter;
    end
end

function [] = validateJitter(x)
    if(ischar(x))
        validatestring(x, {'dynamic'});
    else
        validateattributes(x, {'numeric'}, ...
            {'scalar','nonempty','real','nonnan','nonnegative', '<=', 1})
    end
end

function [] = validateNormCoef(x)
    if(ischar(x))
        validatestring(x, {'auto'});
    else
        validateattributes(x, {'numeric'}, ...
            {'vector','nonempty','real','nonnan','positive'})
    end
end

