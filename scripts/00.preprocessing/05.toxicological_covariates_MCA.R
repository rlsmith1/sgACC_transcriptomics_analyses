

###############################################################################

# Run Multiple Correspondence Analysis on toxicology data across samples
# to reduce dimensionality of drug use data for GRCCA

###############################################################################

# Setup ---------------------------------------------------------------

### SET DIRECTORIES ###

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/supplement/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

### LOAD LIBRARIES ###

library(tidyverse)
library(DESeq2)
library(MASS)

### LOAD DATA ###

## Cleaned covariates
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in 00.clean_covariates.R)


## Define toxicological covariates
drug_covariates <- c(
    "smoker", "nicotine_cotinine", "alcohol", "opioids", "cannabinoids", "beta_blockers",
    "major_stimulants_cocaine_included", "minor_stimulants", "anticholinergics", "antidepressants",
    "anti_epileptics", "anti_histamines", "antipsychotics", "mood_stabilizers", "sedative_hypnotic_anxiolitics", "benzos",
    "other_psychotropic_drug", "non_psychiatric"
)



# Run MCA -----------------------------------------------------------------


## Select drugs for analysis
df_drugs <- df_covariates_numeric %>% 
    ungroup %>% 
    dplyr::select(sample, all_of(drug_covariates)) %>% 
    distinct %>% 
    mutate_all( ~ as.character(.x)) %>% 
    column_to_rownames("sample")

## Run MCA
n_factors <- ncol(df_drugs) # 17 known drugs

set.seed(20240308)
mca.res <- mca(df_drugs %>%
                   # remove NAs (for now)
                   mutate_all(~ifelse(is.na(.x), 0, .x)) %>% 
                   mutate_all(~factor(.x)),
               nf = n_factors
)

## Calculate eigenvalue of each dimension
df_eig <- tibble(
    dim = 1:n_factors,
    var = (mca.res$d^2)*100,
    cum_var = cumsum(var)
)

## Extract variable loadings
df_var_loadings <- mca.res$cs %>% 
    as.data.frame %>% 
    dplyr::rename_all( ~ paste0("dim", .x)) %>% 
    rownames_to_column("variable") %>% 
    as_tibble() %>% 
    separate(variable, into = c("drug", "use"), sep = "\\.") %>% 
    mutate(use = ifelse(use == 1, "yes", "no"))

## Extract individual loadings on each dimension
df_ind_loadings <- mca.res$rs %>%
  as.data.frame %>%
  dplyr::rename_all( ~ paste0("dim", .x)) %>%
  rownames_to_column("sample") %>%
  as_tibble()


# Correlate each covariate with each MCA dimension -----------------------------------------

df_mca_covar_cor <- expand_grid(
    covariate = df_covariates_numeric %>% dplyr::select(-sample) %>% colnames,
    dimension = df_ind_loadings %>% dplyr::select(-sample) %>% colnames
) %>% 
    mutate(
        pearsons_r = map2(
            .x = covariate,
            .y = dimension,
            .f = ~ cor.test(df_covariates_numeric[[.x]], df_ind_loadings[[.y]])$estimate
        ),
        p_val = map2(
            .x = covariate,
            .y = dimension,
            .f = ~ cor.test(df_covariates_numeric[[.x]], df_ind_loadings[[.y]])$p.value
        )
    ) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    mutate(p_adj = p.adjust(p_val, method = "fdr"))


# Save objects ------------------------------------------------------------


## Save for figures
save(df_eig, df_var_loadings, df_mca_covar_cor,
     file = paste0(analysis_objects_dir, "drug_MCA_results.Rdata")
)


## Save for downstream analyses
save(df_var_loadings, df_ind_loadings,
     file = paste0(base_dir, "objects/drug_MCA_results.Rdata")
)


# count NAs for each covariate --------------------------------------------

df_nas <- df_covariates_numeric %>% 
  is.na %>% 
  colSums() %>% 
  enframe %>% 
  dplyr::rename_all(~c("covariate", "n")) %>% 
  mutate(perc_na = n/nrow(df_covariates_numeric)) %>% 
  arrange(-n)

df_nas

# DEPR: covariate PCA ---------------------------------------------------------------------

pca <- df_vsd_regress %>% column_to_rownames("sample") %>% t() %>% 
  prcomp(center = TRUE, scale = TRUE)

df_samples_pca <- pca$x %>% 
  as.data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  as_tibble

df_covariates_pca <- pca$rotation %>% 
  as.data.frame %>% 
  rownames_to_column("sample") %>% 
  as_tibble %>% 
  left_join(df_covariates) 


# PLOT COVARIATES
df_covariates_pca %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = dx), shape = 21) +
  scale_fill_manual(values = dx_colors)


# DEPR: Drug tetrachoric correlation -----------------------------------------------


df_drugs <- df_covariates_numeric %>% 
  ungroup %>% 
  dplyr::select(sample, all_of(drug_covariates)) %>% 
  distinct %>% 
  mutate_all( ~ as.character(.x)) %>% 
  column_to_rownames("sample")

ncol(df_drugs) # 17 drugs

# run tetrachoric correlation
het.res <- hetcor(df_drugs)
het.mat <- het.res$cor
het.mat[is.na(het.mat)] <- 0

# plot drug correlations
drug_order <- het.mat[hclust(dist(het.mat))$order,] %>% rownames

het.mat %>% 
  as.data.frame %>% 
  rownames_to_column("drug1") %>% 
  as_tibble() %>% 
  pivot_longer(2:ncol(.), names_to = "drug2", values_to = "correlation") %>% 
  mutate(drug1 = factor(drug1, levels = drug_order),
         drug2 = factor(drug2, levels = drug_order)
  ) %>% 

  ggplot(aes(x = drug1, y = drug2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
  labs(x = "", y = "", title = "Drug tertachoric correlation matrix") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/covariates/drug_correlation_matrix", .x),
                width = 6, height = 5)
)


# DEPR: Drug factor analysis -----------------------------------------------

# # determine optimal number of factors for factor analysis
# df_factanal <- tibble()
# for (i in 2:10) {
#   
#   n_factors <- i
#   fa.res <- fa(r = cor(het.mat), nfactors = n_factors, rotate = "varimax", SMC = FALSE, fm = "minres")
#   var_explained <- (fa.res$Vaccounted)["Cumulative Var", ] %>% max
#   df_tmp <- tibble(nfactors = n_factors, cum_var = var_explained)
#   df_factanal <- df_factanal %>% bind_rows(df_tmp)
#   
# }
# df_factanal
# 
# # plot nfactors vs cumulative variance explained
# df_factanal %>% 
#   ggplot(aes(x = nfactors, y = cum_var)) +
#   geom_point() +
#   labs(x = "Number of factors", y = "Cumulative variance explained",
#        title = "Selecting number of factors")
# 
# map(
#   .x = c(".png", ".pdf"),
#   .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/covariates/drug_FA_variance_explained", .x),
#                 width = 4, height = 3)
# )
# 
# # select final number of factors
# fa.res <- fa(r = cor(het.mat), nfactors = 5, rotate = "varimax", SMC = FALSE, fm = "minres")
# 
# # plot
# fa.res$score.cor
# fa.res$r.scores
# fa.res$loadings
# fa.res$weights
# 
# loadings <- fa.res$loadings
# df_fa_res <- data.frame(matrix(as.numeric(loadings), attributes(loadings)$dim, dimnames = attributes(loadings)$dimnames)) %>% 
#   rownames_to_column("drug") %>% 
#   as_tibble()
# 
# df_fa_res %>% 
#   pivot_longer(2:ncol(.), names_to = "factor", values_to = "loading") %>% 
#   mutate(factor = str_remove(factor, "MR") %>% as.numeric %>% as.factor) %>% 
#   mutate(drug = factor(drug, levels = drug_order)) %>% 
#   mutate(label = ifelse(abs(loading) > 0.30, round(loading, 3), NA)) %>% 
#   
#   ggplot(aes(x = factor, y = drug, fill = loading)) +
#   geom_tile() +
#   geom_text(aes(label = label), size = 3) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
#   labs(x = "Factor", y = "", 
#        title = "Drug loading on each factor", caption = "Label = abs(loading) > 0.30")
# 
# map(
#   .x = c(".png", ".pdf"),
#   .f = ~ ggsave(paste0(base_dir, "outputs/figures/drug_FA_loadings", .x),
#                 width = 8, height = 5)
# )
# 
# # save for GRCCA
# save(df_fa_res, file = paste0(base_dir, "objects/drug_FA_results.RDS"))

