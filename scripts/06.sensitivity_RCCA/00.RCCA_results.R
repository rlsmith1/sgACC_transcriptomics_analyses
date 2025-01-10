
########################################################################################

# Load gene-level RCCA results

########################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/supplement/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WGCNA_res.R"))


# Load data ---------------------------------------------------------------

## Identify CCA directory
n_mca_dim <- 8
type <- "GENES/"
project_dir <- paste0(prefix,
                      "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "_", 
                      n_mca_dim, "MCA_regressBrainWeight_WITHGRAY/")
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)

## Read in RCCA analysis results
analysis_type <- "rcca"
VARx <- "0.1_1"
l_rcca_res <- f_load_cca_res(
    cca_directory = cca_dir,
    level = "gene",
    analysis_type = analysis_type,
    mu = 0.1,
    VARx = VARx,
    include_Cmat = FALSE
)

## Assign results to tibbles; negate everything so SCZ is positive
df_results_rcca <- l_rcca_res$model_results
df_lvs_rcca <- l_rcca_res$latent_variables %>% mutate_at(vars(c(lvx, lvy)), ~ -.x)
df_y_res_rcca <- l_rcca_res$y_res %>% mutate_at(vars(c(weight, z_score, pearsons_r)), ~ -.x)
df_x_res_rcca <- l_rcca_res$x_res %>% 
    mutate_at(vars(c(weight, z_score, pearsons_r)), ~ -.x) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    left_join(df_modules_filt %>% dplyr::select(-c(mod_set, color)), by = join_by(ensembl_gene_id)) %>% 
    dplyr::select(ensembl_gene_id, gene_symbol, module, everything()) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol))


# Save for plotting -------------------------------------------------------


# Figures
save(df_lvs_rcca, df_y_res_rcca, df_x_res_rcca,
     file = paste0(analysis_objects_dir, "RCCA_results.Rdata")
)

# Downstream analyses
save(df_lvs_rcca, df_y_res_rcca, df_x_res_rcca,
     file = paste0(base_dir, "objects/RCCA_results.Rdata")
)

