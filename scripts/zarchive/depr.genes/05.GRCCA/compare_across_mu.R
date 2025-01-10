
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/full_analysis_scripts/scripts/load_WGCNA_res.R"))
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/05.GRCCA/f_load_CCA_res.R"))

# Load data ---------------------------------------------------------------

## Identify CCA directory
n_mca_dim <- 8
type <- "GENES/"
project_dir <- paste0(prefix,
                      "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "_", 
                      n_mca_dim, "MCA_regressBrainWeight_WITHGRAY/")
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)

## Identify mu values
cca_dir_files <- list.files(paste0(cca_dir, "framework"))
VARx <- "0.1_1"
analysis_string <- paste0("grcca_permutation_VARx", VARx)
cca_dir_files_grcca_VARx <- cca_dir_files[str_detect(cca_dir_files, analysis_string)]
mu_values <- str_remove(cca_dir_files_grcca_VARx, paste0(analysis_string, "_L2xNfeat_groupL2x")) %>%
    str_remove("_noCmat_allEffects")

## Load across mu values
l_grcca_res <- map(
    .x = mu_values,
    .f = ~ f_load_cca_res(
        cca_directory = cca_dir,
        analysis_type = "grcca",
        mu = .x,
        VARx = VARx,
        include_Cmat = FALSE
    )
    
)
names(l_grcca_res) <- mu_values



# Compare XY correlations and variance explained --------------------------

df_model_results <- map(l_grcca_res, 
    ~ .x$model_results %>% 
        mutate(mu = .x$mu, .before = 1)
) %>% 
    bind_rows()

## Plot XY correlation p-values across varX & mu values
df_model_results %>% 
    ggplot(aes(x = varx %>% as.factor(), y = -log10(pval), color = mu)) +
    geom_point() +
    geom_line(aes(group = mu)) +
    geom_hline(aes(yintercept = -log10(0.05)), lty = 2)

## Plot XY correlation across varX & mu values
df_model_results %>% 
    ggplot(aes(x = varx %>% as.factor(), y = correl, color = mu)) +
    geom_point() +
    geom_line(aes(group = mu))



# Compare Y res -----------------------------------------------------------

df_y_res_by_mu <- map(l_grcca_res, 
    ~ .x$y_res %>% 
        mutate(mu = .x$mu, .before = 1)
) %>% 
    bind_rows(.)
   
df_y_res_by_mu %>% 
    ggplot(aes(x = reorder(covariate, pearsons_r), y = pearsons_r, fill = mu)) +
    geom_col(position = "dodge")

l_grcca_res$`0.1` %>% View()


# Compare X res -----------------------------------------------------------


df_x_res_by_mu <- map(l_grcca_res, 
                      ~ .x$x_res %>% 
                          mutate(mu = .x$mu, .before = 1)
) %>% 
    bind_rows(.)

top_genes <- df_x_res_by_mu %>% 
    dplyr::filter(mu == 0.1) %>% 
    top_n(n = 50, wt = abs(pearsons_r)) %>% 
    pull(ensembl_gene_id)

df_x_res_by_mu %>% 
    dplyr::filter(ensembl_gene_id %in% top_genes) %>% 
    ggplot(aes(x = reorder(ensembl_gene_id, pearsons_r), y = pearsons_r, fill = mu)) +
    geom_col(position = "dodge")



