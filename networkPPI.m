function [Wppi] = networkPPI(Sites, STRING_Proteins, STRING_PPI)
    nSiteInitial = height(Sites);

    validSites = Sites.ProteinIndex > 0;
    Sites = Sites(validSites, :);
    
    [~, ib] = ismember(Sites.Protein, STRING_Proteins);
    validSites = Sites.MonoPhosphorylated & (ib > 0);
    Sites = Sites(validSites, :);
    [~, ib] = ismember(Sites.Protein, STRING_Proteins);
    Sites.STRINGproteinIndex = ib;
    
    nSite = height(Sites);
    nSTRINGprotein = length(STRING_Proteins);
    Msite2stringprotein = sparse((1:nSite), Sites.STRINGproteinIndex, ...
            true, nSite, nSTRINGprotein);
    
    Wppi = mapnetwork(Msite2stringprotein, STRING_PPI, Msite2stringprotein');
    Wppi = logical(Wppi - diag(diag(Wppi)));
    
    Mpsite2site = sparse((1:nSite), Sites.SiteIndex, ...
            true, nSite, nSiteInitial);
    Wppi = mapnetwork(Mpsite2site', Wppi, Mpsite2site);
    Wppi = logical(Wppi - diag(diag(Wppi)));
end








