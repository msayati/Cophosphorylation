
function [Wsharedkinase] = networkSharedKinase(Sites)
    nSiteInitial = height(Sites);

    load('in/kinase_substrate_processed.mat');
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
    Mkinase2site = mapnetwork(KSm, Mpspsite2site);

    Wsharedkinase = mapnetwork(Mkinase2site', Mkinase2site);
    Wsharedkinase = logical(Wsharedkinase - diag(diag(Wsharedkinase)));
    
    Mpsite2site = sparse((1:nSite), Sites.SiteIndex, ...
            true, nSite, nSiteInitial);
    Wsharedkinase = mapnetwork(Mpsite2site', Wsharedkinase, Mpsite2site);
    Wsharedkinase = logical(Wsharedkinase - diag(diag(Wsharedkinase)));
end