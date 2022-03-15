function [Wintra_protein] = networkIntraProtein(Sites)
    nSiteInitial = height(Sites);

    validSites = Sites.ProteinIndex > 0;
    Sites = Sites(validSites, :);
    
    nSite = height(Sites);
    nProtein = max(Sites.ProteinIndex);
    Msite2protein = sparse((1:nSite), Sites.ProteinIndex, ...
            true, nSite, nProtein);
    Wintra_protein = mapnetwork(Msite2protein, Msite2protein');
    Wintra_protein = logical(Wintra_protein - diag(diag(Wintra_protein)));
    
    Mpsite2site = sparse((1:nSite), Sites.SiteIndex, ...
            true, nSite, nSiteInitial);
    Wintra_protein = mapnetwork(Mpsite2site', Wintra_protein, Mpsite2site);
    Wintra_protein = logical(Wintra_protein - diag(diag(Wintra_protein)));
end