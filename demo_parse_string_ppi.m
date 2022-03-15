Start = tic();
%%
load('in/string_ppi_network.mat');
%%
PPI_high_conf = PPI > 700;

filePath = 'in/string_ensp_to_gene_name.tab';
ds = tabularTextDatastore(filePath, 'FileExtensions', '.tab');
ds.TextscanFormats = repmat({'%q'}, 1, length(ds.VariableNames));
Mapping = ds.readall();

Mapping.PrimaryGene = regexprep(Mapping.Gene, ' [^.]*', '');
Mapping.PrimaryEnsp = regexprep(Mapping.ENSP, ',[^.]*', '');

[STRING_Proteins, ~, prot_index] = unique(Mapping.PrimaryGene);
[ENSP, ~, ensp_index] = unique(Mapping.PrimaryEnsp);

Mensp2stringp = logical(sparse(ensp_index, prot_index, ...
        1, length(ENSP), length(STRING_Proteins)));
    
[b, ib] = ismember(ENSP, Proteins);
if(nnz(~b) > 0)
   error('Invalid ensp identifier.'); 
end
STRING_PPI = PPI_high_conf(ib, ib);
STRING_PPI = STRING_PPI | speye(size(STRING_PPI));
STRING_PPI = mapnetwork(Mensp2stringp', STRING_PPI, Mensp2stringp);
STRING_PPI = logical(STRING_PPI - diag(diag(STRING_PPI)));

save('out/string_network.mat', 'STRING_Proteins', 'STRING_PPI');







