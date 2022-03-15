function [Wptmcode] = networkPTMsigdb(Sites, ...
                        PTMprotein, PTMpsite, Mptmsite2ptmsite)
    nSiteInitial = height(Sites);

    nPTMprotein = length(PTMprotein);
    nPTMpsite = height(PTMpsite);
    nPosition = max(PTMpsite.Position);

    [~, ib] = ismember(Sites.Protein, PTMprotein);
    validSites = Sites.MonoPhosphorylated & (ib > 0) ...
            & (Sites.Position <= nPosition);
    Sites = Sites(validSites, :);

    [~, ib] = ismember(Sites.Protein, PTMprotein);
    Sites.PTMproteinIndex = ib;
    Sites.PSIndex = sub2ind([nPTMprotein nPosition], ...
            Sites.PTMproteinIndex, Sites.Position);

    [~, ib] = ismember(Sites.PSIndex, PTMpsite.PSIndex);
    Sites.PTMpsiteIndex = ib;
    validSites = (ib > 0);
    Sites = Sites(validSites, :);

    nSite = height(Sites);
    Msite2ptmsite = sparse((1:nSite), Sites.PTMpsiteIndex, ...
            true, nSite, nPTMpsite);

    Wptmcode = mapnetwork(Msite2ptmsite, Mptmsite2ptmsite, Msite2ptmsite');
    Wptmcode = logical(Wptmcode - diag(diag(Wptmcode)));
    
    Mpsite2site = sparse((1:nSite), Sites.SiteIndex, ...
            true, nSite, nSiteInitial);
    Wptmcode = mapnetwork(Mpsite2site', Wptmcode, Mpsite2site);
    Wptmcode = logical(Wptmcode - diag(diag(Wptmcode)));
end