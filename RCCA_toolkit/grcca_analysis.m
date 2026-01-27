function grcca_analysis(mu)

%----- Analysis

% Set path for analysis
set_path;

% Project folder
base_dir = '/data/smithral/sgacc_wgcna/RCCA_toolkit/';
%type = 'OVERALLsubs/';
%project_dir = '08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC_CVq1_sft2_minSize35_cutHeight0.988_8MCA_regressDim2_WITHGray';
type = 'GENES/'
project_dir = '08Mar2024_GENES_qSVAgeSexRaceGC_sft3_minSize40_cutHeight0.98_8MCA_regressBrainWeight_WITHGRAY'
cfg.dir.project = strcat(base_dir, type, project_dir);

% Machine settings
cfg.machine.name = 'grcca';
cfg.machine.param.L2y = 0; % lambda y
cfg.machine.param.groupL2y = 0; % mu y
cfg.machine.svd.varx = 0.99; % 1;
cfg.machine.svd.vary = 1;
cfg.machine.param.type = 'factorial';

% --> maximize variance explained
cfg.machine.param.name = {'VARx' 'L2x' 'groupL2x' 'L2y' 'groupL2y'};
cfg.machine.param.VARx = 0.1:0.1:1;

% --> set x hyperparameters
X = load(strcat(base_dir, type, project_dir, '/data/X.mat'));
n_features = size(X.X, 2);
lambda = 1 - (1 / n_features); % feature-level regularization is dependent on number of features
cfg.machine.param.L2x = lambda; % lambda x
cfg.machine.param.groupL2x = mu; % mu x (on a quasi-logarithmic scale, so 0.1 is set to almost 0 & 0.9999 is almost 1)

% Data settings
cfg.data.conf = 1; % turn off (to 0) if not regressing confounders (no C.mat)

% Framework settings
cfg.frwork.name = 'permutation'; % optimize variance explained
%cfg.frwork.name = 'holdout'; % optimize hyperparameters
cfg.frwork.split.nout = numel(cfg.machine.param.VARx); % use when optimizing variance explained
%cfg.frwork.split.nout = 1;
%cfg.frwork.nlevel = 1; % turn off to get all associative effects
cfg.frwork.flag = strcat('_VARx0.1_1_L2xNfeat_groupL2x', num2str(mu), '_noCmat_allEffects_TEST'); % use this to differentiate between analyses with same name

% Environment settings
cfg.env.comp = 'cluster'; % cluster or local
%cfg.env.save.svd = 1; % turn off to save disc space
%cfg.env.save.tableHeading = {'set' 'correl' 'pval' 'varx' 'npcax' 'l2x' 'groupl2x'};

% Statistical inference settings
cfg.stat.nperm = 1000;
cfg.stat.nboot = 1000;

% Update cfg with defaults
cfg = cfg_defaults(cfg);

% Run analysis
main(cfg);

% Clean up analysis files to save disc space
cleanup_files(cfg);

end