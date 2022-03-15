rng(1, 'twister'); % For reproducibility
%%
% datasetname = 'AD';
datasets = {'BC', 'BC-PDX', 'CPTAC', 'OV', ...
    'OV-PDX', 'CRC', 'LC', 'AD', 'RPE'};
dataset_labels = {'BC1', 'BC2', 'BC3', 'OC1', ...
    'OC2', 'CRC', 'LC', 'AD', 'RPE'};
% datasets = {'OV', 'OV-PDX'};
logtransform = true;
scatter = false;
measure = 'corr';
permstyle = 'nodes';
network = 'intraprotein';
network2 = 'permute';
% network = 'ptmcode-inter';
% network2 = 'permute';
proximity = 20;
proximity2 = 20;
motif_length = 2;
motif_length2 = 1;
filtersites = true;
nBin = 'auto';
% nBin = 100;
xlimits = 'auto';
% xlimits = [0 3];
xlimits = [-1 1];
% legendpos = 'best';
legendpos = 'southeast';
figure_no = 2;
%%
experiments = {'intra', 'inter'};
nExperiment = numel(experiments);
if(nExperiment ~= 2)
    error('Number of experiments should be 2.');
end
    
nDataset = numel(datasets);
D_all = cell(nExperiment, nDataset);
Dp_all = cell(nExperiment, nDataset);
for iDataset = 1:nDataset
    datasetname = datasets{iDataset};
    load(sprintf('in/%s.mat', datasetname));
    if(logtransform)
        X = log2(X);
    end
    nSite = size(X, 1);
    nSample = size(X, 2);

    Wintra_protein = networkIntraProtein(Sites);
    networkColorIndex = 1;
    switch(network)
        case 'all'
            W = true(nSite, nSite);
            networkname = 'AllPairs';
        case 'intraprotein'
            W = Wintra_protein;
            networkname = 'IntraProtein';
        case 'proximity'
            W = networkSequenceProximity(Sites, Wintra_protein, proximity);
            networkname = sprintf('SeqProximity%d', proximity);
        case 'sharedkinase'
            W = networkSharedKinase(Sites);
            networkname = 'SharedKinase';
            networkColorIndex = 1;
        case 'sharedkinase-intra'
            W = networkSharedKinase(Sites);
            W = W & Wintra_protein;
            networkname = 'SharedKinase-Intra';
        case 'sharedkinase-inter'
            W = networkSharedKinase(Sites);
            W = logical(W - (W & Wintra_protein));
            networkname = 'SharedKinase-Inter';
        case 'motif'
            W = networkMotif(Sites, motif_length);
            networkname = sprintf('Motif%d', motif_length);
        case 'motif-intra'
            W = networkMotif(Sites, motif_length);
            W = W & Wintra_protein;
            networkname = sprintf('Motif%d-Intra', motif_length);
        case 'motif-inter'
            W = networkMotif(Sites, motif_length);
            W = logical(W - (W & Wintra_protein));
            networkname = sprintf('Motif%d-Inter', motif_length);
        case 'ppi'
            load('in/string_ppi_network.mat');
            W = networkPPI(Sites, STRING_Proteins, STRING_PPI);
            networkname = 'PPI';
            networkColorIndex = 2;
        case 'ptmcode'
            load('in/human_phosphosite_network_combined.mat');
            W = networkPTMcode(Sites, PTMprotein, PTMpsite, Mptmsite2ptmsite);
            networkname = 'PTMcode';
            networkColorIndex = 3;
        case 'ptmsigdb'
            load('out/ptmsigdb_network.mat');
            W = networkPTMsigdb(Sites, PTMprotein, PTMpsite, Mptmsite2ptmsite_ptmsigdb);
            networkname = 'PTMsigdb';
            networkColorIndex = 4;
        case 'ptmcode-intra'
            load('in/human_phosphosite_network_combined.mat');
            W = networkPTMcode(Sites, PTMprotein, PTMpsite, Mptmsite2ptmsite);
            W = W & Wintra_protein;
            networkname = 'PTMcode-Intra';
        case 'ptmcode-inter'
            load('in/human_phosphosite_network_combined.mat');
            W = networkPTMcode(Sites, PTMprotein, PTMpsite, Mptmsite2ptmsite);
            W = logical(W - (W & Wintra_protein));
            networkname = 'PTMcode-Inter';
        otherwise
            error('Invalid network.');
    end
    W_base = W;
    X_base = X;
    Sites_base = Sites;

    for iExperiment = 1:nExperiment
        X = X_base;
        Sites = Sites_base;
        experiment = experiments{iExperiment};
        switch(lower(experiment))
            case 'intra'
                W = W_base & Wintra_protein;
                network2 = 'intraprotein';
            case 'inter'
                W = logical(W_base - (W_base & Wintra_protein));
                network2 = 'permute';
            otherwise
                error('Invalid experiment name.');
        end
        switch(network2)
            case 'permute'
                Wp = W_base;
            case 'all'
                Wp = true(nSite, nSite);
                networkname2 = 'AllPairs';
            case 'intraprotein'
                Wp = Wintra_protein;
                networkname2 = 'IntraProtein';
            case 'proximity'
                Wp = networkSequenceProximity(Sites, Wintra_protein, proximity2);
                networkname2 = sprintf('SeqProximity%d', proximity2);
            case 'sharedkinase'
                Wp = networkSharedKinase(Sites);
                networkname2 = 'SharedKinase';
            case 'motif'
                Wp = networkMotif(Sites, motif_length);
                networkname2 = sprintf('Motif%d', motif_length2);
            case 'motif-intra'
                Wp = networkMotif(Sites, motif_length);
                Wp = Wp & Wintra_protein;
                networkname2 = sprintf('Motif%d-Intra', motif_length2);
            case 'motif-inter'
                Wp = networkMotif(Sites, motif_length);
                Wp = logical(Wp - (Wp & Wintra_protein));
                networkname2 = sprintf('Motif%d-Inter', motif_length2);
            case 'ptmcode'
                load('in/human_phosphosite_network_combined.mat');
                Wp = networkPTMcode(Sites, PTMprotein, PTMpsite, Mptmsite2ptmsite);
                networkname2 = 'PTMcode';
            case 'ptmcode-intra'
                load('in/human_phosphosite_network_combined.mat');
                Wp = networkPTMcode(Sites, PTMprotein, PTMpsite, Mptmsite2ptmsite);
                Wp = Wp & Wintra_protein;
                networkname2 = 'PTMcode-Intra';
            case 'ptmcode-inter'
                load('in/human_phosphosite_network_combined.mat');
                Wp = networkPTMcode(Sites, PTMprotein, PTMpsite, Mptmsite2ptmsite);
                Wp = logical(Wp - (Wp & Wintra_protein));
                networkname2 = 'PTMcode-Inter';
            otherwise
                error('Invalid network.');
        end

        measure_ = measure;
        if(strcmpi(measure, 'absdifpos'))
            X = Sites.Position;
            validsites = ~isnan(Sites.Position);
            X = X(validsites, :);
            W = W(validsites, validsites);
            Wp = Wp(validsites, validsites);
            Sites = Sites(validsites, :);
            nSite = size(X, 1);
            measure_ = 'absdif';
        end

        if(filtersites)
            validsites = (sum(W, 1) > 0) | (sum(Wp, 1) > 0);
            X = X(validsites, :);
            W = W(validsites, validsites);
            Wp = Wp(validsites, validsites);
            Sites = Sites(validsites, :);
            nSite = size(X, 1);
        end

        if(strcmpi(network2, 'permute'))
            switch(permstyle)
                case 'samples'
                    Xp = zeros(size(X));
                    for iSite = 1:nSite
                        Xp(iSite, :) = X(iSite, randperm(nSample));
                    end
                case 'sites'
                    Xp = zeros(size(X));
                    for iSample = 1:nSample
                        Xp(:, iSample) = X(randperm(nSite), iSample);
                    end
                case 'nodes'
                    perm = randperm(nSite);
                    Wp = Wp(perm, perm);
                    Xp = X;
%                     Xp = X(randperm(nSite), :);
                otherwise
                    error('Invalid permutation style.');
            end
            networkname2 = sprintf('perm%s', permstyle);
            networkname2_lgd = sprintf('%s - Permuted (%s)', networkname, permstyle);
            networkname_lgd = sprintf('%s - Original', networkname);
        else
            Xp = X;
            networkname_lgd = networkname;
            networkname2_lgd = networkname2;
        end

        D = assessCophos(X, W, measure_);
        Dp = assessCophos(Xp, Wp, measure_);
        D_all{iExperiment, iDataset} = D';
        Dp_all{iExperiment, iDataset} = Dp';
    end
end
% D = horzcat(D_all{:});
% Dp = horzcat(Dp_all{:});
%%
switch(measure)
    case 'absdif'
        ylabeltext = 'Mean Phosphorylation Difference';
    case 'absdifpos'
        ylabeltext = 'Mean Position Difference';
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
description = 'Co-phosphorylation of sites';

violinOpts = struct();
violinOpts.ViolinAlpha = 0.6;
violinOpts.NormCoef = 2.1;
% violinOpts.GroupNames = datasets;
violinOpts.Boundaries = [-1.1 1.1];
violinOpts.Clipped = false;
violinOpts.ShowMean = false;
% violinOpts.ViolinPosition = 'left';
violinOpts.Width = 1;
violinOpts.Jitter = 'dynamic';
violinOpts.MeanLineStyle = '-';
violinOpts.MeanLineWidth = 3;
violinOpts.MeanLine = 0.9;
% violinOpts.LineWidth = 1.25;
% violinOpts.JitterPosition = 'left';
% violinOpts.Scatter = true;

Info = table();
Info.Dataset = dataset_labels';
Info.MeanDif = zeros(nDataset, 1);
Info.MeanDifCI = zeros(nDataset, 2);
Info.TStat = zeros(nDataset, 1);
Info.PValues = zeros(nDataset, 1);
Info.KSpval = zeros(nDataset, 1);
for iDataset = 1:nDataset
    % Left - Right
   Dleft = Dp_all{1, iDataset};
   Dright = Dp_all{2, iDataset};
   [h, p, ci, stats] = ttest2(Dleft, Dright);
   mdif = mean(Dleft) - mean(Dright);
   Info.MeanDif(iDataset) = mdif;
   Info.MeanDifCI(iDataset, :) = ci;
   Info.TStat(iDataset) = stats.tstat;
   Info.PValues(iDataset) = p;
   Info.SampleNleft(iDataset) = length(Dleft);
   Info.SampleNright(iDataset) = length(Dright);
   [h, p, ksstat] = kstest2(Dleft, Dright);
   Info.KSpval(iDataset) = p;
end
figure_path = sprintf('out/figures/paper_output/figure2/');
if(~exist(figure_path, 'dir')); mkdir(figure_path); end
save([figure_path, 'figure2a_stats_data.mat'], 'Info');



figure(figure_no);
clf();
hold('on');
colors = get(gca, 'ColorOrder');
for iExperiment = 1:2
    if(iExperiment == 1)
        violinOpts.ViolinPosition = 'left';
        violinOpts.JitterPosition = 'left';
%         violinOpts.GroupColors = [0.598, 0.58, 0.58];
        violinOpts.GroupColors = [0.925 0.67 0.55];
    else
        scatter = false;
        violinOpts.ViolinPosition = 'right';
        violinOpts.JitterPosition = 'right';
        violinOpts.GroupColors = [0.3 0.3 0.9];
    end
    for iDataset = 1:nDataset

%         permColor = [0.925 0.67 0.55];
%         permColor = brighten([0.598, 0.58, 0.58], 0.1);
        permColor = [0.7 0.7 0.7];
%         colorMain = [0.2 0.2 0.9];
        colorMain = colors(networkColorIndex, :);
%         h1 = violinplot(D_all(iExperiment, iDataset), 'XShift', iDataset-1, ...
%             'GroupColors', colorMain, 'Scatter', false, ...
%             'MarkerColor', [0.05 0.05 0.9], violinOpts, ...
%             'Front', 'scatter', 'ViolinAlpha', 0.6, ...
%             'MeanLineColor', 'groupcolor');
        h2 = violinplot(Dp_all(iExperiment, iDataset), 'XShift', iDataset-1, ...
            'GroupColors', permColor, 'Scatter', false, violinOpts, ...
            'LineWidth', 0.1, 'ViolinAlpha', 0.5, 'ShowMean', true, ...
            'Darkening', 0.6);
%         violinplot(D_all(iExperiment, iDataset), 'XShift', iDataset-1, ...
%                     'GroupColors', colorMain, 'Scatter', false, ...
%                     'MarkerColor', [0.05 0.05 0.9], violinOpts, ...
%                     'Violin', false, 'MeanLineColor', 'darken', ...
%                     'Darkening', 0.15, 'ShowMean', true);
    end
end
% h1 = violinplot(D_all, ...
%     'GroupColors', [0.2 0.2 0.9], violinOpts);
% h2 = violinplot(Dp_all, ...
%     'GroupColors', [0.925 0.67 0.55], violinOpts);
hold('off');
set(gca, 'XTick', (1:nDataset) - 0.5, 'XTickLabel', dataset_labels);
% ylabel(sprintf('%s (%s)', ylabeltext, fctext));
ylabel('Co-phosphorylation');
% title(description);
% legend([h1(1) h2(1)], {sprintf('%s', networkname_lgd), ...
%     sprintf('%s', networkname2_lgd)}, 'Location', legendpos);
set(gca, 'FontSize', 13);
set(gcf, 'Color', [1 1 1]);
xlim('auto');
ylim('auto');
pause(0.02);
xlim([0 nDataset]);
ylim([-1 1]);

set(gcf, 'Position', [0, 0, 1280, 400]);
movegui('center');

figure_path = sprintf('out/figures/paper_output/figure2/');
figure_name = sprintf('violinplot_intra_inter_%s_%s', fctext, measure);
if(~exist(figure_path, 'dir')); mkdir(figure_path); end
export_fig([figure_path, figure_name, '.png'], '-m3');
%%
return;
%%
for iExperiment = 1:2
    D = zeros(1, nDataset);
    N = zeros(1, nDataset);
    for iDataset = 1:nDataset
       m1 = mean(D_all{iExperiment, iDataset});
       m2 = mean(Dp_all{iExperiment, iDataset});
       n1 = length(D_all{iExperiment, iDataset});
       d = m1 - m2;
       D(iDataset) = d;
       N(iDataset) = n1;
    end
    if(iExperiment == 1) % Intra-Protein
        intraD = mean(D);
        intraN = log10(mean(N));
    else                 % Inter-Protein
        interD = mean(D);
        interN = log10(mean(N));
    end
end

if(isinf(intraN))
    intraN = 0;
    intraD = 0;
end

P = [intraD intraN interD interN]
labels = {{'Intra-Protein', 'Co-Phosphorylation'}, {'Intra-NumEdges'}, ...
    {'Inter-Protein', 'Co-Phosphorylation'}, {'Inter-NumEdges'}};
%     limits = [0 0 0 0; 1 1 0.3 0.3];
% limits = [0 0 0 0; 1 1 0.3 0.3];
limits = [0 0 0 0; 0.3 6 0.3 6];
colorP = colors(networkColorIndex, :);

reordering = [1 4 3 2];

P = P(reordering);
labels = labels(reordering);
limits = limits(:, reordering);

figure(3);
clf();
spider_plot(P, ...
    'AxesLabels', labels, ...
    'AxesLimits', limits, ...
    'AxesInterval', 3, ...
    'Color', colorP, ...
    'FillOption', 'on', ...
    'AxesPrecision', 1, ...
    'AxesFontSize', 15, ...
    'Log10Labeling', [false, true, false, true], ...
    'AxesLabelsOffset', 0.12, ...
    'HorzCenterAll', true, ...
    'LabelFontSize', 15);
set(gcf, 'Color',  [1 1 1]);
set(gcf, 'Position', [0 0 680 680]);
movegui('center');

figure_path = sprintf('out/figures/dual_violinplot/');
figure_name = sprintf('%s_spider_intra_inter_%s_%s', ...
    lower(networkname), fctext, measure);
if(~exist(figure_path, 'dir')); mkdir(figure_path); end
export_fig([figure_path, figure_name, '.png'], '-m3');


%%

