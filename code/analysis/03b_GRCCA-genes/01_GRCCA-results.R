
#==============================================================================#
# Characterize gene-level GRCCA results (gene weights & enrichments)
#==============================================================================#

## Load GRCCA results (output from toolkit script grcca_analysis.m on cluster)
## Run permutation analysis to test the significance of rx association in 
## each diagnostic group


## Setup
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(scripts_dir, "load_WGCNA_res.R"))


# Load data ---------------------------------------------------------------

## Identify CCA directory
dx_groups_to_analyze <- c("Control", "BD", "MDD", "SCZ") # select diagnostic groups of interest
n_mca_dim <- 8
type <- "GENES/"
project_dir <- paste0(prefix,
                      "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "_", 
                      paste0(dx_groups_to_analyze, collapse = ""), "_",
                      n_mca_dim, "MCA_regressBrainWeight_WITHGRAY/")
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)


## Read in GRCCA analysis results (f_load_cca_res defined in functions directory)
analysis_type <- "grcca"
VARx <- "0.1_1"
l_grcca_res <- f_load_cca_res(
    cca_directory = cca_dir,
    level = "gene",
    analysis_type = analysis_type,
    mu = 0.1,
    VARx = VARx,
    include_Cmat = FALSE,
    rename_covariates = TRUE
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



# Test the significance of the structure correlation alignment with expression in each diagnostic group ------------------

## Combine structure correlation and expression data
df_rx_expr <- df_vsd_regress %>% 
    pivot_longer(contains("ENSG"), names_to = "ensembl_gene_id", values_to = "expression_value") %>% 
    mutate(dx = case_when(
        str_detect(sample, "control") ~ "Control",
        str_detect(sample, "bipolar") ~ "BD",
        str_detect(sample, "mdd") ~ "MDD",
        str_detect(sample, "schizo") ~ "SCZ"
    ) %>% factor(levels = names(dx_colors))) %>% 
    group_by(ensembl_gene_id, dx) %>% 
    summarise(mean_expr = mean(expression_value)) %>% 
    left_join(df_x_res)


## Calculate the actual (observed) correlation between rx and mean expression in each diagnostic group
df_rx_expr_cor <- df_rx_expr %>% 
    group_by(dx) %>% 
    group_modify( ~ {
        
        test <- cor.test(.x$mean_expr, .x$pearsons_r)
        
        tibble(
            correlation = test$estimate,
            p_value = test$p.value
        )
        
    }) %>% 
    ungroup() %>% 
    mutate(p_adj = p.adjust(p_value, method = "BH"))

## Permute the rx vector and re-calculate the correlation
n_perm <- 10000
df_rx_expr_cor_perm <- tibble()
set.seed(0206)
for (i in 1:n_perm) {
    
    if (i %% 100 == 0) {print(paste0("Running permutation ", i))}
    
    # 1: Shuffle structure correlations within module
    df_perm <- df_rx_expr %>%
        group_by(module) %>%
        mutate(permuted_r = sample(pearsons_r, replace = FALSE)) %>%
        ungroup()
    
    # For each diagnostic group, correlate permuted_r with mean_expr
    df_tmp <- df_perm %>% 
        group_by(dx) %>% 
        group_modify( ~ {
            
            test <- cor.test(.x$mean_expr, .x$permuted_r)
            
            tibble(
                correlation = test$estimate,
                p_value = test$p.value
            )
            
        }) %>% 
        ungroup() %>% 
        mutate(
            p_adj = p.adjust(p_value, method = "BH"),
            perm = i
        )
    
        # Append to tibble
        df_rx_expr_cor_perm <- df_rx_expr_cor_perm %>% bind_rows(df_tmp)
    
}

## Calculate permutation p-values
df_rx_expr_cor_perm_pvals <- df_rx_expr_cor_perm %>% 
    left_join(df_rx_expr_cor %>% select(dx, obs_corr = correlation), by = "dx") %>%
    group_by(dx) %>%
    summarise(
        perm_p = (sum(abs(correlation) >= abs(obs_corr)) + 1) / (n_perm + 1)
    ) %>% 
    left_join(df_rx_expr_cor)


# Save -------------------------------------------------------


save(df_lvs, df_results, df_y_res, df_x_res, df_rx_expr_cor_perm_pvals,
     file = paste0(objects_dir, "GRCCA_results.Rdata")
)


