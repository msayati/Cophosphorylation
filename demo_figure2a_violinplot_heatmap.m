inFolder = 'out/figures/paper_output/figure2/';
load([inFolder, 'figure2a_stats_data.mat']);
%%
M = Info.MeanDif';
T = Info.TStat;
CI = Info.MeanDifCI';

C = min(abs(CI(1, :)), abs(CI(2, :)));
signChange = sign(CI(1, :)) ~= sign(CI(2, :));
C(signChange) = 0;

[nRow, nColumn] = size(C);

labels = cell(size(C));
for iRow = 1:nRow
    for iColumn = 1:nColumn
        mVal = M(iRow, iColumn);
        mText = sprintf('\\fontsize{%d} %.2f', 11, mVal);
        ciText = sprintf('\\fontsize{%d} [%.2f, %.2f]', 10, CI(1, iColumn), CI(2, iColumn));
        
        label = {mText, ciText};
        labels{iRow, iColumn} = label;
    end
end

colorMin = [1 1 1];
colorMax = [0.98 0.5 0.1];
% colorMax = [0.99 0.57 0.45];
colorMax = brighten(colorMax, -0.1);


aconst = 0.05;
fq = @(a) 1 - 2./(1+2.^(a/aconst));
fcolor = @(q) color_spacing_continuous(q, [0 1], [1 1 1; colorMax]);
fcomb = @(a) fcolor(fq(a));

colors = fcomb(C(:));
colors = reshape(colors, nRow, nColumn, 3);

figure(1);
clf();
% cbbox(C);
cbbox(C, 'Labels', labels, ...
    'cMin', 0, 'cMax', 0.5, ...
    'colorMin', colorMin, 'colorMax', colorMax, ...
    'boxWidthRatio', 0.9, ...
    'Colors', colors);
% set(gca, 'FontSize', 12);
% set(gca, 'XTickLabel', Info.Dataset);

set(gca, 'Visible', 'off');

set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position', [0 0 1000 60]);
movegui('center');
export_fig([inFolder, 'figure2a_stats_bar.png'], '-m2');
