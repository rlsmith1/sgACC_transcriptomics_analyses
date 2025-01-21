% set dir
base_dir = '/data/smithral/sgacc_wgcna/RCCA_toolkit/';
%type= 'GENES/';
type = 'OVERALLsubs/';
project_dir = '08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC_CVq1_sft2_minSize35_cutHeight0.988_ControlBDMDDSCZ_8MCA_regressDim2_broadRareSCZBDMDDASD_sigGRCCAfdr0_05';

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

% convert GroupsX
group = importdata(strcat(base_dir, type, project_dir, '/data/GroupsX.txt'));
group = group.data;
save(strcat(base_dir, type, project_dir, '/data/GroupsX.mat'), 'group');
