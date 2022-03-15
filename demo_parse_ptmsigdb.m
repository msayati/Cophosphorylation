Start = tic();
%%
load('in/data_PTMsigDB_all_sites_v1.9.0_human.mat');
%%
T = dataPTMsigDB;
validRows = ismember(T.category, {'PATH-NP', 'PATH-WP'});
T = T(validRows, :);
a = cellstr(split(T.siteuniprot, ';'));
T.uniprotacc = a(:, 1);
T.site = a(:, 2);
[T.residue, T.position] = splitSiteResidue(T.site);
a = cellstr(split(T.siteannotation, '_'));
T.protein = a(:, 1);

T.identifier = strcat(T.protein, '_', T.site);

[pathways, ~, pathway_index] = unique(T.signature);
T.pathway_index = pathway_index;
nPathway = length(pathways);
T.direction_up = strcmpi(T.sitedirection, 'u');
T.direction_down = strcmpi(T.sitedirection, 'd');

% T = T(T.direction_up, :);

[PTMprotein, ~, protein_index] = unique(T.protein);
T.protein_index = protein_index;

[sites, table_index, site_index] = unique(T.identifier);
T.site_index = site_index;

nSite = length(sites);
nPathway = length(pathways);
nProtein = length(PTMprotein);
nPosition = max(T.position);

PTMpsite = table();
PTMpsite.Identifier = sites;
PTMpsite.Protein = T.protein(table_index);
PTMpsite.ProteinIndex  = T.protein_index(table_index);
PTMpsite.Site  = T.site(table_index);
PTMpsite.Position  = T.position(table_index);
PTMpsite.PSIndex = sub2ind([nProtein nPosition], ...
            PTMpsite.ProteinIndex, PTMpsite.Position);

Mpathway2ptmsite = sparse(T.pathway_index, T.site_index, ...
        true, nPathway, nSite);
        
Mptmsite2ptmsite_ptmsigdb = mapnetwork(Mpathway2ptmsite', Mpathway2ptmsite);
Mptmsite2ptmsite_ptmsigdb = Mptmsite2ptmsite_ptmsigdb | speye(size(Mptmsite2ptmsite_ptmsigdb));
        
save('out/ptmsigdb_network.mat', 'PTMprotein', 'PTMpsite', 'Mptmsite2ptmsite_ptmsigdb');







