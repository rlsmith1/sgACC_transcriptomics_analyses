% set dir
base_dir = '/Users/smithral/Documents/PhD/projects/sgacc_wgcna_grcca/RCCA_toolkit/';
type = 'GENES/';
project_dir = '08Mar2024_GENES_qSVAgeSexRaceGC_sft3_minSize40_cutHeight0.98_8MCA_regressBrainWeight_WITHGRAY';

% convert X
X = importdata(strcat(base_dir, type, project_dir, '/data/X.txt'));
X = X.data;
save(strcat(base_dir, type, project_dir, '/data/X.mat'), 'X');

% convert Y
Y = importdata(strcat(base_dir, type, project_dir, '/data/Y.txt'));
Y = Y.data;
save(strcat(base_dir, type, project_dir, '/data/Y.mat'), 'Y');

% convert C
C = importdata(strcat(base_dir, type, project_dir, '/data/C.txt'));
C = C.data;
save(strcat(base_dir, type, project_dir, '/data/C.mat'), 'C');

% convert XGroup
group = importdata(strcat(base_dir, type, project_dir, '/data/XGroup.txt'));
group = group.data;
save(strcat(base_dir, type, project_dir, '/data/XGroup.mat'), 'group');
