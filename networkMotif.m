function [Wmotif] = networkMotif(Sites, d)
    nSiteInitial = height(Sites);
    load('in/phosphosites_processed.mat');
    
    [proteins] = unique(Sites.Protein);
    validSites = Sites.ProteinIndex > 0;
    Sites = Sites(validSites, :);

    nSite = height(Sites);
    nProtein = length(proteins);
    nPosition = max(max(Sites.Position), max(S.Position));

    Msite2protein = sparse((1:nSite), Sites.ProteinIndex, ...
        true, nSite, nProtein);
    
    [~, ib] = ismember(S.Gene, proteins);
    Mpspsite2protein = sparse(1:height(S), ib + 1, ...
            true, height(S), nProtein + 1);
    Mpspsite2protein = Mpspsite2protein(:, 2:end);
    Mpspsite2site1 = mapnetwork(Mpspsite2protein, Msite2protein');

    Mpspsite2position = sparse(1:height(S), S.Position, ...
            true, height(S), nPosition);

    pos = Sites.Position;
    pos(isnan(pos)) = 0;
    Msite2position = sparse(1:nSite, pos + 1, ...
            true, nSite, nPosition + 1);
    Msite2position = Msite2position(:, 2:end);

    Mpspsite2site2 = mapnetwork(Mpspsite2position, Msite2position');
    Mpspsite2site = Mpspsite2site1 & Mpspsite2site2;
    
    
    validSites = sum(Mpspsite2site, 1) == 1;
    Mpspsite2site(:, ~validSites) = 0;
    [pspindices, siteindices] = ind2sub(size(Mpspsite2site), find(Mpspsite2site));

    Sites.PSPIndex = zeros(nSite, 1);
    Sites.PSPIndex(siteindices) = pspindices;
    Sites.Sequence = cell(nSite, 1);
    Sites.Sequence(siteindices) = S.Sequence(pspindices);
    Sites.Binding = cell(nSite, 1);
    Sites.Binding(:) = {''};
    Sites.Binding(siteindices) = cellfun(@(x) upper(x(8-d:8+d)), ...
        S.Sequence(pspindices), 'UniformOutput', false);

    [bindingsites, ~, ib] = unique(Sites.Binding);
    Msite2binding = sparse(1:nSite, ib, ...
            true, nSite, length(bindingsites));
    Msite2binding(:, 1) = 0;
    Wmotif = mapnetwork(Msite2binding, Msite2binding');
    Wmotif = logical(Wmotif - diag(diag(Wmotif)));
    
    Mpsite2site = sparse((1:nSite), Sites.SiteIndex, ...
            true, nSite, nSiteInitial);
    Wmotif = mapnetwork(Mpsite2site', Wmotif, Mpsite2site);
    Wmotif = logical(Wmotif - diag(diag(Wmotif)));
end