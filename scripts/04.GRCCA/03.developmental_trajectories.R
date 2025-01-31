
################################################################################

# Characterize normative developmental expression of GRCCA results deciles

################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/Fig4_characterizeGRCCA/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/oad_WGCNA_res.R"))

## Load GRCCA data
load(paste0(base_dir, "objects/GRCCA_results.Rdata")) # df_lvs, df_y_res, df_x_res



# Developmental trajectories (data) ----------------------------------------------

load(paste0(base_dir, "objects/psychencode_development_expr_data.Rdata")) # df_psychencode_metadata, df_psychencode_expr

## Identify samples from region of interest
df_sample_window <- df_psychencode_metadata %>% 
    dplyr::filter(region_broad == "frontal_cortex")
frontal_cortex_samples <- df_sample_window%>% 
    pull(id)

## Filter expression data for samples and genes of interest
df_psychencode_expr_filt <- df_psychencode_expr %>% 
    dplyr::filter(GENE %in% ensembl_gene_ids) %>% 
    dplyr::select(GENE, all_of(frontal_cortex_samples)) %>% 
    dplyr::rename("ensembl_gene_id" = "GENE")

## Combine developmental expression data with GRCCA results (split into deciles)
df_psychencode_grcca <- df_psychencode_expr_filt %>% 
    pivot_longer(starts_with("HSB"), names_to = "id", values_to = "dev_expression") %>% 
    left_join(df_sample_window, by = join_by(id)) %>% 
    dplyr::select(window, id, regioncode, ensembl_gene_id, dev_expression) %>% 
    left_join(df_x_res %>% 
                  mutate(decile = ntile(-pearsons_r, 10) %>% factor)
    )

## Save for figures
save(df_psychencode_grcca,
     file = paste0(analysis_objects_dir, "GRCCA_developmental_expression.RDS")
)

