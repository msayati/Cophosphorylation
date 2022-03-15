function [] = cbbox(C, varargin)
    p = inputParser;
    validX = @(x) validateattributes(x, {'cell', 'numeric', 'logical'}, ...
        {'2d', 'nonempty', 'real'});
    validScalar = @(x) validateattributes(x, {'logical', 'numeric'}, ...
        {'scalar','nonempty','real','nonnan'});
    validColorArray = @(x) validateattributes(x, {'numeric'}, ...
        {'2d', 'nonempty', 'real', 'nonnan', 'ncols', 3});
    addRequired(p, 'C', validX);
    addParameter(p, 'Labels', [], @iscell);
    addParameter(p, 'cMin', min(C(:)), validScalar);
    addParameter(p, 'cMax', max(C(:)), validScalar);
    addParameter(p, 'colorMin', [1 1 1], validColorArray);
    addParameter(p, 'colorMax', [1 0 0], validColorArray);
    addParameter(p, 'boxWidthRatio', 1, validScalar);
    addParameter(p, 'boxHeightRatio', 1, validScalar);
    addParameter(p, 'Colors', [], @isnumeric);
    parse(p, C, varargin{:});
    param = p.Results;
    checkUsingDefaults = @(p,varname) any(strcmp(p.UsingDefaults,varname));
    
    [nRow, nColumn] = size(C);
    
    if(checkUsingDefaults(p, 'Labels'))
        labels = cell(size(C));
        for iRow = 1:nRow
            for iColumn = 1:nColumn
                cValue = C(iRow, iColumn);
                labels{iRow, iColumn} = num2str(cValue, 2);
            end
        end
        param.Labels = labels;
    end
   
    if(checkUsingDefaults(p, 'Colors'))
        colorMin = param.colorMin;
        colorMax = param.colorMax;
        cols = color_spacing_continuous(C(:), [param.cMin; param.cMax], [colorMin; colorMax]);
        param.Colors = reshape(cols, nRow, nColumn, 3);
    end
    
    box_width_ratio = param.boxWidthRatio;
    box_height_ratio = param.boxHeightRatio;
    
    row_height = 1 ./ nRow;
    col_width = 1./ nColumn;
    
    xcPos = zeros(nColumn, 1);
    ycPos = zeros(nRow, 1);
    
    txtOptions = struct();
    txtOptions.HorizontalAlignment = 'center';
    txtOptions.VerticalAlignment = 'middle';
    txtOptions.FontSize = 12;
    for iRow = 1:nRow
        for iColumn = 1:nColumn
            row_gap = (1 - box_height_ratio) * 0.5;
            col_gap = (1 - box_width_ratio) * 0.5;
            x = 0 + (col_gap + iColumn - 1) * col_width;
            y = 0 + (row_gap + iRow - 1) * row_height;
            w = col_width * (1 - col_gap);
            h = row_height * (1 - row_gap);
            xcenter = x + w/2;
            ycenter = y + h/2;
            xcPos(iColumn) = xcenter;
            ycPos(iRow) = ycenter;
            pos = [x y w h];
            color = param.Colors(iRow, iColumn, :);
            rectangle('Position',pos, 'FaceColor', color)
            
            txt = param.Labels{iRow, iColumn};
            text(xcenter, ycenter, txt, txtOptions);
            
        end
    end
    set(gca, 'XTick', xcPos, 'XTickLabel', []);
    set(gca, 'YTick', ycPos, 'YTickLabel', []);
    
    
    
    
end

