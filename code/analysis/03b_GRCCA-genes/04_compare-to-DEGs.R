
#==============================================================================#
# Determine the relationship between DEG test statistics and GRCCA structure correlations
#==============================================================================#

## Characterize the relationship between DGE L2FC and GRCCA rx
## using permutation testing


## Setup
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(scripts_dir, "load_WGCNA_res.R"))

load(paste0(objects_dir, "GRCCA_results.Rdata")) # df_lvs, df_y_res, df_x_res


## Combine GRCCA and DEG results
df_de_grcca <- df_x_res %>% 
    left_join(df_limma_res, by = join_by(ensembl_gene_id, gene_symbol))


## Calculate the observed correlation between rx and DEG L2FC
cor.test(df_de_grcca$pearsons_r, df_de_grcca$log_fc)


## Run permutation testing (resampling within module) to determine the significance of the correlation
n_perm <- 10000
df_limma_grcca_cor_perm <- tibble()
set.seed(0406)
for (i in 1:n_perm) {
    
    if (i %% 100 == 0) {print(paste0("Running permutation ", i))}
    
    # Permute LFC within module
    df_tmp <- df_de_grcca %>% 
        group_by(module) %>%
        mutate(permuted_r = sample(pearsons_r, replace = FALSE))
    
    # Calculate correlation between permuted data and Nirmala's LFC
    test <- cor.test(df_tmp$permuted_r, df_tmp$log_fc)
    correlation <- test$estimate
    p_value <- test$p.value
    
    # Append to tibble
    df_limma_grcca_cor_perm <- df_limma_grcca_cor_perm %>% 
        bind_rows(
            tibble(
                perm = i,
                correlation = correlation,
                p_value = p_value
            )
        )
    
}


## Calculate permutation p-values
de_limma_grcca_perm_p <- df_limma_grcca_cor_perm %>% 
    mutate(obs_corr = cor.test(df_de_grcca$pearsons_r, df_de_grcca$log_fc)$estimate) %>%
    summarise(
        perm_p = (sum(abs(correlation) >= abs(obs_corr)) + 1) / (n_perm + 1)
    )
# Pperm = 1*10^-4


