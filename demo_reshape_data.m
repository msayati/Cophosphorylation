%% BC
load('in/data0/BC_data.mat');
X = BC_phospho;
Sites = parseSiteInfo(BC_sites);
save('in/BC.mat', 'X', 'Sites');
%% BC-PDX
load('in/data0/BC_PDX_phosphodata.mat');
[p, s] = splitSiteName(phosphoproteins, '-');
s2 = regexprep(s, '[a-z]', ',');
ps = cellstr(join([p, s2], '-'));
X = phosphodata;
Sites = parseSiteInfo(ps); % To be fixed
save('in/BC-PDX.mat', 'X', 'Sites');
%% CP-TAC
load('in/data0/CPTAC_data.mat');
% [p, s] = splitSiteName(genes, '-');
% s2 = regexprep(s, '[a-z]', ',');
% ps = cellstr(join([p, s2], '-'));
X = 2.^ phosphodata;
Sites = parseSiteInfo(genes);
save('in/CPTAC.mat', 'X', 'Sites');
%% OV
load('in/data0/OV_data.mat');
[p, s] = splitSiteName(OV_sites', '-');
s2 = regexprep(s, '[a-z]', ',');
ps = cellstr(join([p, s2], '-'));
X = OV_phospho;
Sites = parseSiteInfo(ps);
save('in/OV.mat', 'X', 'Sites');
%% OV-PDX
load('in/data0/OV_PDX_phosphodata.mat');
X = phosphodata;
s2 = regexprep(sites', '[a-z]', ',');
ps = cellstr(join([genes', s2], '-'));
Sites = parseSiteInfo(ps);
save('in/OV-PDX.mat', 'X', 'Sites');
%% CRC
load('in/data0/CRC_data_unique.mat');
valids = sum(isinf(phosphodata), 2) == 0;
X = phosphodata(valids,:);
ps1 = phosphoproteins(valids);
[p, s] = splitSiteName(ps1, '-');
s2 = regexprep(s, '([A-Za-z])(\d)', ',$1$2');
s3 = regexprep(s2, '^.', '');
ps = cellstr(join([p, s3], '-'));
Sites = parseSiteInfo(ps);
save('in/CRC.mat', 'X', 'Sites');
%% LC
load('in/data0/LungCancer_data.mat');
X = phosphodata;
Sites = parseSiteInfo(phosphosites); 
save('in/LC.mat', 'X', 'Sites');
%% AD
load('in/data0/AlzheimerData.mat');
[p, s] = splitSiteName(AD_namaes, '-', false);
s2 = regexprep(s, ';', ',');
s3 = regexprep(s2, '-cleavage.*', ' cleavage');
ps = cellstr(join([p, s3], '-'));
X = AD;
Sites = parseSiteInfo(ps);
save('in/AD.mat', 'X', 'Sites');
%% RPE
load('in/data0/RPE_data.mat');
X = phosphodata;
Sites = parseSiteInfo(phosphosites);
save('in/RPE.mat', 'X', 'Sites');

