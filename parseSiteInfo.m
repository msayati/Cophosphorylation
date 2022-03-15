function [ Sites ] = parseSiteInfo( sites )
    sites = reshape(sites, [], 1);
    nSite = length(sites);
    
    [protein, site] = splitSiteName(sites, '-');
    protein(strcmpi(protein, 'NA')) = {''};
    [u, ~, protein_indices] = unique(protein);
    if(isempty(u{1}))
        protein_indices(protein_indices==1) = 0;
    end
    
    [~, n] = splitSites(site, ',');
    
    monoPSites = n == 1;
    [residue, position] = splitSiteResidue(site(monoPSites), ',');
    
    Sites = table();
    Sites.SiteIndex = (1:nSite)';
    Sites.Protein = protein;
    Sites.ProteinIndex = protein_indices;
    Sites.Site = site;
    Sites.NumSites = n;
    Sites.MonoPhosphorylated = n == 1;
    Sites.Residue = cell(nSite, 1);
    Sites.Residue(monoPSites) = residue;
    Sites.Position = nan(nSite, 1);
    Sites.Position(monoPSites) = position;
end

function [ site, n ] = splitSites( A, delimiter )
    [site, n] = cellfun(@(x) funSplit(x, delimiter), A, ...
        'UniformOutput', false);
    n = cell2mat(n);
end

function [site, n] = funSplit(x, delimiter)
    c = strsplit(x, delimiter);
    n = length(c);
    site = c{1};
end

function [ residue, position ] = splitSiteResidue( A, delimiter )
    [residue, position] = cellfun(@funSplitResidue, A, ...
        'UniformOutput', false);
    position = cell2mat(position);
end

function [residue, position] = funSplitResidue(x)
    residue = upper(x(1));
%     position = x(regexp(x, '\d'));
    position = str2double(x(regexp(x, '\d')));
end


