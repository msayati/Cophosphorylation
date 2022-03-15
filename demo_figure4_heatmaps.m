inFolder = 'out/figures/paper_output/figure4/';
outFolder = inFolder;
isLeft = 1;

%%
networks = {'sharedkinase', 'ppi', 'ptmcode', 'ptmsigdb'};
network_labels = {'SharedKinase', 'PPI', 'Co-Evolution', 'Pathways'};
nDataset = 9;
nNetwork = length(networks);

M = zeros(nDataset, nNetwork);
T = zeros(nDataset, nNetwork);
CImin = zeros(nDataset, nNetwork);
CImax = zeros(nDataset, nNetwork);
C = zeros(nDataset, nNetwork);

for iNetwork = 1:nNetwork 
    network = networks{iNetwork};
    load([inFolder, 'figure4_', network, '_stats_data.mat']);
    if(isLeft)
        Infox = InfoLR.('Left');
        lrText = 'left';
    else
        Infox = InfoLR.('Right');
        lrText = 'right';
    end
    
    M(:, iNetwork) = Infox.MeanDif;
    T(:, iNetwork) = Infox.TStat;
%     T = Info.TStat;

    
    ciMin = Infox.MeanDifCI(:, 1);
    ciMax = Infox.MeanDifCI(:, 2);
    CImin(:, iNetwork) = ciMin;
    CImax(:, iNetwork) = ciMax;
    
    c = min(abs(ciMin), abs(ciMax));
    signChange = sign(ciMin) ~= sign(ciMax);
    signNegative = ciMax < 0;
    c(signChange) = 0;
    c(signNegative) = c(signNegative) * -1;
    C(:, iNetwork) = c;
end

%%
[nRow, nColumn] = size(C);

labels = cell(size(C));
for iRow = 1:nRow
    for iColumn = 1:nColumn
        mVal = M(iRow, iColumn);
        ciMin = CImin(iRow, iColumn);
        ciMax = CImax(iRow, iColumn);
        mText = sprintf('\\fontsize{%f} %.2f', 15, mVal);
        ciText = sprintf('\\fontsize{%f} [%.2f, %.2f]', 13, ciMin, ciMax);
        
        label = {mText, ciText};
        if(isnan(mVal)); label = 'NA'; end
        labels{iRow, iColumn} = label;
    end
end

colorMinus = [0 0 1];
colorZero = [1 1 1];
colorMax = [0.98 0.1 0.1];

defaultcolors = getdefaultcolors();

% aconst = 0.1;
% fq = @(a) 1 - 2./(1+2.^(a/aconst));
% fcolor = @(q) color_spacing_continuous(q, [-1 0 1], [colorMinus; colorZero; colorMax]);
% fcomb = @(a) fcolor(fq(a));
% colors = fcomb(C(:));
% colors = reshape(colors, nRow, nColumn, 3);

colors = zeros(nRow, nColumn, 3);
for iNetwork = 1:nNetwork
    colorPlus = defaultcolors(iNetwork, :);
    colorPlus = brighten(colorPlus, 0.5);
    aconst = 0.05;
    fq = @(a) 1 - 2./(1+2.^(a/aconst));
    fcolor = @(q) color_spacing_continuous(q, [0 1], [colorZero; colorPlus]);
    
%     fcolor = @(q) color_spacing_continuous(q, [-1 0 1], [colorMinus; colorZero; colorMax]);
    fcomb = @(a) fcolor(fq(a));
    colors(:, iNetwork, :) = fcomb(abs(C(:, iNetwork)));
end

xticklabels = network_labels;
for iNetwork = 1:nNetwork
    color = defaultcolors(iNetwork, :);
    color = color * (1 - 0);
    color = brighten(color, -0.25);
    L = relative_luminance(color);
    if((L >= 0.35))
        color = color * (1 - L*0.3);
    end
    label = network_labels{iNetwork};
    xticklabels{iNetwork} = sprintf('\\color[rgb]{%f, %f, %f}\\bf\\fontsize{15}{%s}', color, label);
    
end


figure(1);
clf();
% cbbox(C);
cbbox(C, 'Labels', labels, ...
    'boxWidthRatio', 0.9, ...
    'boxHeightRatio', 1, ...
    'Colors', colors);
set(gca, 'FontSize', 15);
set(gca, 'XTickLabel', xticklabels);

set(gca, 'YDir', 'reverse');

% set(gca, 'YTickLabel', Info.Dataset);

% set(gca, 'Visible', 'off');

set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position', [0 0 600 600]);
movegui('center');
export_fig([outFolder, 'figure4_stats_bar_', lrText, '.png'], '-m2');
