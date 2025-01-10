
## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/updated_figures/figures/Fig4_characterizeGRCCA/")
tables_dir <- paste0(base_dir, "outputs/updated_figures/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/updated_figures/objects/Fig4_characterizeGRCCA/")

source(paste0(base_dir, "scripts/full_analysis_scripts/scripts/load_WGCNA_res.R"))

## Load GRCCA data
load(paste0(base_dir, "objects/GRCCA_results.Rdata")) # df_lvs, df_y_res, df_x_res

## Load DEGs
load(paste0(base_dir, "objects/DE_results.RDS")) # df_de_res

## Combine DE & GRCCA res
df_grcca_de_res <- df_de_res %>% 
    left_join(df_x_res)



# Astrocytes? -------------------------------------------------------------

cell_types[["Astro"]]
df_astrocyte_genes <- df_grcca_de_res %>% 
    dplyr::filter(ensembl_gene_id %in% cell_types[["Astro"]])

df_astrocyte_genes %>% 
    mutate(label = ifelse(abs(log2fold_change) > 0.075, gene_symbol, "") ) %>% 
    ggplot(aes(x = log2fold_change, y = pearsons_r)) +
    geom_point() +
    geom_text_repel(aes(label = label)) +
    stat_cor()

high_FC_low_r <- df_astrocyte_genes %>% 
    dplyr::filter(abs(log2fold_change) > 0.075) %>% 
    pull(ensembl_gene_id)
df_vsd_regress %>% 
    dplyr::select(sample, all_of(high_FC_low_r)) %>% 
    left_join(df_covariates) %>% 
    pivot_longer(contains("ENSG"), names_to = "ensembl_gene_id", values_to = "expression_value") %>% 
    left_join(df_ensembl_to_symbol, by = join_by(ensembl_gene_id)) %>% 

    ggplot(aes(x = mood_stabilizers, y = expression_value)) +
    geom_point(position = position_jitter(width = 0.2), size = 0.75) +
    geom_boxplot(fill = "transparent", outlier.shape = NA) +
    facet_wrap(vars(gene_symbol))


df_x_res %>% dplyr::filter(gene_symbol == "NRXN1")

df_vsd_regress %>% 
    dplyr::select(sample, all_of("ENSG00000179915")) %>% 
    left_join(df_covariates) %>% 
    pivot_longer(contains("ENSG"), names_to = "ensembl_gene_id", values_to = "expression_value") %>% 
    left_join(df_ensembl_to_symbol, by = join_by(ensembl_gene_id)) %>% 
    
    ggplot(aes(x = anti_histamines, y = expression_value)) +
    geom_point(position = position_jitter(width = 0.2), size = 0.75) +
    geom_boxplot(fill = "transparent", outlier.shape = NA) +
    facet_wrap(vars(gene_symbol))



# SNAP paper loadings -----------------------------------------------------


df_ling_factors <- read_xlsx("~/Downloads/41586_2024_7109_MOESM6_ESM.xlsx", sheet = "Gene loadings")

df_lf4 <- df_ling_factors %>% 
    dplyr::select(id, PEER_4) %>% 
    separate(id, into = c("gene_symbol", "cell_type"), sep = "_")

df_lf4_grcca_de <- df_lf4 %>% 
    left_join(df_grcca_de_res, by = join_by(gene_symbol))


df_lf4_grcca_de %>% 
    dplyr::filter(ensembl_gene_id %in% c(cell_types[["Neuro-Ex"]], cell_types[["Neuro-In"]], cell_types[["Astro"]])) %>% 
    mutate(cell_type = ifelse(cell_type %in% c("gabaergic", "glutamatergic"), "neuron", cell_type)) %>% 
    
    ggplot(aes(x = PEER_4, y = pearsons_r)) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_cor() +
    facet_wrap(vars(cell_type), scales = "free")




# snRNAseq ----------------------------------------------------------------


df_ex_de <- read_xlsx("~/Downloads/inline-supplementary-material-13.xlsx") %>% 
    clean_names


df_ex_de %>% 
    left_join(df_grcca_de_res, by = join_by(gene_symbol)) %>% 
    dplyr::filter(!is.na(ensembl_gene_id)) %>% 
    
    ggplot(aes(x = adjusted_p_value, y = pearsons_r)) +
    geom_point()



