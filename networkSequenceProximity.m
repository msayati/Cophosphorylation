function [Wproximity] = networkSequenceProximity(Sites, Wintra_protein, proximity)
    nSiteInitial = height(Sites);

    validSites = Sites.MonoPhosphorylated;
    Sites = Sites(validSites, :);
    Wintra_protein = Wintra_protein(validSites, validSites);

    nSite = height(Sites);
    nPosition = max(Sites.Position);

    Msite2position = sparse((1:nSite), Sites.Position, ...
            true, nSite, nPosition);
    Mposition2position = sparse(generatePositionMatrix(nPosition, proximity));

    Mseq_dist = mapnetwork(Msite2position, Mposition2position, Msite2position');
    Wproximity = Mseq_dist & Wintra_protein;
    
    Mpsite2site = sparse((1:nSite), Sites.SiteIndex, ...
            true, nSite, nSiteInitial);
    Wproximity = mapnetwork(Mpsite2site', Wproximity, Mpsite2site);
    Wproximity = logical(Wproximity - diag(diag(Wproximity)));
end
