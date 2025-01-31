
########################################################################################

# Compare genes in SCZ WGCNA modules with risk gene lists

########################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/supplement/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WGCNA_res.R"))

### Module eigengene analysis results ###
load(paste0(base_dir, "objects/kME_analysis_res.RDS")) # df_lm_dx
sig_modules <- df_lm_dx %>% 
    dplyr::filter(p_value < 0.05 & term == "SCZ") %>% 
    pull(module)


### Run hypergeometric test of each WGCNA module with each risk gene list ###
modules <- df_modules_filt %>% pull(module) %>% unique
df_mods_risk_gene_hypergeometric <- expand_grid(
    benchmark_list = names(benchmarking_lists),
    module = modules
) %>% 
    mutate(
        hypergeometric_res = map2(
            .x = benchmark_list,
            .y = module,
            .f = ~ f_module_hypergeometric(
                gene_list = benchmarking_lists[[.x]], 
                gene_list_name = .x, 
                wgcna_module = .y, 
                gene_universe = ensembl_gene_ids, 
                module_data = df_modules_filt
            ) %>% 
                dplyr::select(-module)
        )
    ) %>% 
    unnest(cols = c(hypergeometric_res)) %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr"))

## Save for supplemental figure
save(
    df_mods_risk_gene_hypergeometric,
    file = paste0(analysis_objects_dir, "module_risk_gene_hypergeometric.RDS")
)
