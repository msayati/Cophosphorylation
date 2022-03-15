rng(1, 'twister'); % For reproducibility
%%
datasets = {'BC', 'BC-PDX', 'CPTAC', 'OV', ...
    'OV-PDX', 'CRC', 'LC', 'AD', 'RPE'};
dataset_labels = {'BC1', 'BC2', 'BC3', 'OC1', ...
    'OC2', 'CRC', 'LC', 'AD', 'RPE'};
measure = 'corr';
figure(1);
clf();
for iDataset = 1:length(datasets)
    datasetname = datasets{iDataset};
    datasetlabel = dataset_labels{iDataset};
    fprintf('Dataset %d: %s\n', iDataset, datasetname);
    
    subplot(3, 3, iDataset);
    drawPlot(datasetname, datasetlabel, iDataset, measure);
end

set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position',  [0 0 1100, 720]);
movegui('center');
outFolder = 'out/figures/paper_output/figure2/';
if(~exist(outFolder, 'dir')); mkdir(outFolder); end
figure_name = sprintf('combined_%s_errorbar', measure);
export_fig([outFolder, figure_name, '.png'], '-m2');
%%
function [] = drawPlot(datasetname, datasetlabel, iSubplot, measure)
    [iColumn, iRow] = ind2sub([3 3], iSubplot);
    
    logtransform = true;
    alpha = 0.05;
    xlims = [0 3];

    load(sprintf('in/%s.mat', datasetname));
    if(logtransform)
        X = log2(X);
    end

    validSites = Sites.ProteinIndex > 0 & Sites.MonoPhosphorylated;
    Sites = Sites(validSites, :);
    X = X(validSites, :);

    nSite = size(X, 1);
    nSample = size(X, 2);
    nProtein = max(Sites.ProteinIndex);
    nPosition = max(Sites.Position);

    Msite2protein = sparse((1:nSite), Sites.ProteinIndex, ...
            true, nSite, nProtein);
    Wintra_protein = mapnetwork(Msite2protein, Msite2protein');
    Wintra_protein = logical(Wintra_protein - diag(diag(Wintra_protein)));

    % [D] = assessPairwiseMean(X, Wintra_protein);
    D = assess(X, Wintra_protein, measure);
    P = assessPairwiseMean(Sites.Position, Wintra_protein);

    % % Xp = X(randperm(nSite), :);
    % [Dp] = assessPairwiseMean(X, Wintra_protein);
    switch(measure)
        case 'absdif'
            ylabeltext = 'Mean Phosphorylation Difference';
        case 'abscorr'
            ylabeltext = 'Absolute Correlation';
        case 'corr'
            ylabeltext = 'Correlation';
        case 'bicorr'
            ylabeltext = 'Biweight Midcorrelation';
        otherwise
            error('Invalid co-phosphorylation measure.');
    end

    if(logtransform)
        fctext = 'log2-FC';
    else
        fctext = 'FC';
    end
    description = 'Co-phosphorylation of close sites';

    nPoints = 10000;
    xx = zeros(1, nPoints);
    M = zeros(1, nPoints);
    S = zeros(1, nPoints);
    SE = zeros(1, nPoints);
    CONF = zeros(1, nPoints);
    for i = 1:nPoints
        vals = D(P<=i);
        n = numel(vals);
        xx(i) = i;
        M(i) = mean(vals);
        S(i) = std(vals);
        SE(i) = std(vals) ./ sqrt(n);
        CONF(i) = SE(i) * tinv(1 - alpha/2, n-1);
    end

    hold('on');
    shadedErrorBar(xx, M, CONF);
    plot(0, 0);
    hold('off');
    grid();
    if(strcmpi(measure, 'corr'))
        ylim([-1 1]);
    end
    if(strcmpi(measure, 'abscorr'))
        ylim([0 1]);
    end
%     ylabel(sprintf('%s (%s)', ylabeltext, fctext));
%     ylabel(sprintf('%s', ylabeltext, fctext));
    ylabel('Co-phosphorylation');
%     if(iColumn == 1)
%         ylabel('Co-phosphorylation');
%     end
%     xlabel('Sequence Proximity (Number of amino acids)');
%     xlabel('Sequence Proximity');
      xlabel('Distance on sequence');
%     if(iRow == 3)
%         xlabel('Sequence Proximity');
%     end
    title(sprintf('%s', datasetlabel));
%     title(sprintf('%s, Dataset:%s', description, datasetname));
    set(gca, 'XScale', 'log');
    set(gca, 'FontSize', 10);
end

function [D] = assess(X, W, measure)
    A = cophospho(X, 'measure', measure);
    m = size(W, 1);
    [i, j] = ind2sub(size(W), find(triu(W, 1)));
    D = A((i-1).*(m-i/2)+j-i);    
%     D = A(find(triu(W, 1)));
end

