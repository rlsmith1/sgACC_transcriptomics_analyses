
########################################################################################

# Load transcript-level GRCCA results

########################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/Fig5_GRCCAtranscriptsRes/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WTCNA_res.R"))


# Load data ---------------------------------------------------------------

# 08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC_CVq1_sft2_minSize35_cutHeight0.988_ControlBDMDDSCZ_8MCA_ControlBDMDDSCZ_regressBrainWeight_WITHGRAY
# 08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC_CVq1_sft2_minSize35_cutHeight0.988_ControlBDMDDSCZ_8MCA_regressDim2_WITHGRAY


## Identify CCA directory
dx_groups_to_analyze <- c("Control", "BD", "MDD", "SCZ") # select diagnostic groups of interest
n_mca_dim <- 8
type <- "OVERALLsubs/"
project_dir <- paste0(prefix, # date & preprocessing info
                      "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "_", # WTCNA parameters
                      paste0(dx_groups_to_analyze, collapse = ""), "_", # diagnostic groups included
                      n_mca_dim, "MCA_regressDim2_", # covariates included
                      "broadRareSCZBDMDDASD_sigGRCCAfdr0_05/" # define subset of transcripts included (based on gene-level results)
                      #"WITHGRAY" # includes all transcripts
)
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)

## Read in GRCCA analysis results
analysis_type <- "grcca"
VARx <- "0.1_1"
l_grcca_res <- f_load_cca_res(
    cca_directory = cca_dir,
    level = "transcript",
    analysis_type = analysis_type,
    mu = 0.1,
    lambda = "Nfeat",
    VARx = VARx,
    include_Cmat = TRUE,
    rename_covariates = FALSE
)

## Assign results to tibbles; negate everything so SCZ is positive
df_results_transcripts <- l_grcca_res$model_results
df_lvs_transcripts <- l_grcca_res$latent_variables %>% mutate_at(vars(c(lvx, lvy)), ~ -.x)
df_y_res_transcripts <- l_grcca_res$y_res %>% mutate_at(vars(c(weight, z_score, pearsons_r)), ~ -.x)
df_x_res_transcripts <- l_grcca_res$x_res %>% 
    mutate_at(vars(c(weight, z_score, pearsons_r)), ~ -.x) %>% 
    left_join(df_transcript_to_gene, by = join_by(ensembl_transcript_id)) %>% 
    left_join(df_modules_filt %>% dplyr::select(-c(mod_set, color)), by = join_by(ensembl_transcript_id)) %>% 
    dplyr::select(ensembl_transcript_id, transcript_symbol, module, everything()) %>% 
    mutate(transcript_symbol = ifelse(is.na(transcript_symbol), ensembl_transcript_id, transcript_symbol)) %>% 
    distinct()


# Save for plotting -------------------------------------------------------


# Figures
save(df_results_transcripts, df_lvs_transcripts, df_y_res_transcripts, df_x_res_transcripts,
     file = paste0(analysis_objects_dir, "GRCCA_transcript_results.Rdata")
)

# Downstream analyses
save(df_results_transcripts, df_lvs_transcripts, df_y_res_transcripts, df_x_res_transcripts,
     file = paste0(base_dir, "objects/GRCCA_transcript_results.Rdata")
)


