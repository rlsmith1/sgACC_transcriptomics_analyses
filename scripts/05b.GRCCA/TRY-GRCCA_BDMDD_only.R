
################################################################################

# Remove SCZ from Y mat to see if we get a BD or MDD signal?

################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/updated_figures/figures/Fig3_GRCCAres/")
tables_dir <- paste0(base_dir, "outputs/updated_figures/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/updated_figures/objects/Fig3_GRCCAres/")

source(paste0(base_dir, "scripts/full_analysis_scripts/scripts/load_WGCNA_res.R"))


# Load data ---------------------------------------------------------------

## Identify CCA directory
dx_groups_to_analyze <- c("Control", "BD", "MDD") # select diagnostic groups of interest
n_mca_dim <- 8
type <- "GENES/"
project_dir <- paste0(prefix,
                      "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "_", 
                      paste0(dx_groups_to_analyze, collapse = ""), "_",
                      n_mca_dim, "MCA_regressBrainWeight_WITHGRAY/")
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)

## Read in GRCCA analysis results
analysis_type <- "grcca"
VARx <- "0.1_1"
l_grcca_res <- f_load_cca_res(
    cca_directory = cca_dir,
    level = "gene",
    analysis_type = analysis_type,
    mu = 0.1,
    VARx = VARx,
    include_Cmat = FALSE,
    rename_covariates = FALSE
)

## Assign results to tibbles; negate everything so SCZ is positive
df_results <- l_grcca_res$model_results
df_lvs <- l_grcca_res$latent_variables %>% mutate_at(vars(c(lvx, lvy)), ~ -.x)
df_y_res <- l_grcca_res$y_res %>% mutate_at(vars(c(weight, z_score, pearsons_r)), ~ -.x)
df_x_res <- l_grcca_res$x_res %>% 
    mutate_at(vars(c(weight, z_score, pearsons_r)), ~ -.x) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    left_join(df_modules_filt %>% dplyr::select(-c(mod_set, color)), by = join_by(ensembl_gene_id)) %>% 
    dplyr::select(ensembl_gene_id, gene_symbol, module, everything()) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol))


