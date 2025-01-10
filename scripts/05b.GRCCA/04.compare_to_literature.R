
################################################################################

# Relate GRCCA results vector to SCZ literature

################################################################################


### SETUP ###

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/supplement/")
tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WGCNA_res.R"))


## Load GRCCA data
load(paste0(base_dir, "objects/GRCCA_results.Rdata")) # df_lvs, df_y_res, df_x_res


### SNAP PAPER DATA ###

## Load SNAP paper LF loadings
df_loading_factors <- read_xlsx("~/Downloads/41586_2024_7109_MOESM6_ESM.xlsx", sheet = "Gene loadings")

## Subset for LF 4
df_lf4 <- df_loading_factors %>% 
    dplyr::select(id, PEER_4) %>% 
    separate(id, into = c("gene_symbol", "cell_type"), sep = "_")

## Combine with GRCCA results
df_lf4_grcca <- df_lf4 %>% 
    left_join(df_x_res, by = join_by(gene_symbol))


### SAVE FOR FUTURE ANALYSES ###
save(df_lf4_grcca, file = paste0(analysis_objects_dir, "SNAP_LF4.RDS"))



### PLOT ###
df_lf4 %>% 
    left_join(df_de_res, by = join_by(gene_symbol)) %>% 
    left_join(df_lake_cell_type) %>% 
    dplyr::filter((type == "Astro" & cell_type == "astrocyte") |
                      (type == "Neuro-Ex" & cell_type == "glutamatergic") |
                      (type == "Neuro-In" & cell_type == "gabaergic")
    ) %>% 
    
    ggplot(aes(x = PEER_4, y = stat)) +
    geom_vline(aes(xintercept = 0), color = "gray") +
    geom_hline(aes(yintercept = 0), color = "gray") +
    geom_point(aes(color = stat), size = 0.5) +
    geom_smooth(color = "black", method = "lm", linewidth = 0.75) +
    stat_cor(label.sep = "\n", size = 3, 
             label.y.npc = "top", label.x.npc = "right", hjust = 1, vjust = 1) +
    facet_wrap(vars(cell_type), scales = "free") +
    scale_color_gradientn(colors = gene_weight_color_scale, guide = "none") +
    scale_x_continuous(breaks = c(-2000, 0, 2000)) +
    labs(x = "Ling et al LF4", y = "GRCCA structure correlation (r)",
         title = "S12 | GRCCA correlation with SNAP program by cell type")

