
########################################################################

# Correlate gene expression vectors with module eigengene to calculate
# gene association with each module

########################################################################

### SETUP ###

## Set directories
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/WGCNA_res/")

## Load data and functions
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/load_WGCNA_res.R"))

# module eigengene analysis results
load(paste0(base_dir, "objects/kME_analysis_res.RDS"))
sig_modules <- df_lm_dx %>% 
    dplyr::filter(p_value < 0.05 & term == "SCZ") %>% 
    pull(module)

gene_universe <- df_vsd_regress %>% dplyr::select(-sample) %>% colnames
modules <- df_kme %>% pull(module) %>% unique

### RUN LOOP ###

df_gene_module_cor <- tibble()
doParallel::registerDoParallel()
for (m in 1:length(modules)) {
    
    current_module <- modules[m]
    print(paste0("Calculating gene-module correlations for ", current_module))
    
    for (g in 1:length(gene_universe)) {
        
        current_gene <- gene_universe[g]
        
        if (g %% 100 == 0) {print(paste0("gene ", g, " of ", length(gene_universe)))}
        
        df_kme_tmp <- df_kme %>% 
            ungroup %>% 
            filter(module == current_module) %>% 
            dplyr::select(sample, kme)
        
        df_gene_tmp <- df_vsd_regress %>% 
            dplyr::select(sample, all_of(current_gene)) %>% 
            dplyr::rename_at(2, ~ "gene")
        
        df_kme_gene_tmp <- df_kme_tmp %>% 
            left_join(df_gene_tmp, by = join_by(sample))
        
        cor_tmp <- cor.test(df_kme_gene_tmp$kme, df_kme_gene_tmp$gene)$estimate
        p_tmp <- cor.test(df_kme_gene_tmp$kme, df_kme_gene_tmp$gene)$p.value
        
        df_tmp <- tibble(
            module = current_module,
            ensembl_gene_id = current_gene,
            pearsons_r = cor_tmp,
            p_value = p_tmp
        )
        
        df_gene_module_cor <- df_gene_module_cor %>% bind_rows(df_tmp)
        
    }
    
}

df_gene_module_cor <- df_gene_module_cor %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr"))

save(df_gene_module_cor, file = paste0(base_dir, "objects/gene_expr_module_kME_correlation.RDS"))
