
#==============================================================================#
# Run WGCNA pipeline 
## (here on gene-level data; same applies to transcript-level)
#==============================================================================#

## WGCNA pipeline:
## (1) Select soft-thresholding power
## (2) Calculate adjacency matrix & TOM matrix --> cluster
## (3) Assign modules

#----Recommended to run on HPC---#


# Setup -------------------------------------------------------------------

## Set directories and load data

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/load_genes.R"))


# (1) Select soft thresholding power ------------------------------------------


## Convert normalized & regressed counts to transposed matrix
m_vsd_regress <- df_vsd_regress %>% column_to_rownames("sample")


## Choose a set of soft-thresholding powers to test
powers <- 1:16


## Call the network topology analysis function
sft <- pickSoftThreshold(m_vsd_regress, powerVector = powers, verbose = 5)


## Convert results to a tibble
df_sft <- sft$fitIndices %>% as_tibble %>% clean_names

## Save
save(df_sft, file = paste0(base_dir, "objects/", prefix, "_SFT.RDS")) # for downstream analyses
save(df_sft, file = paste0(analysis_objects_dir, prefix, "_SFT.RDS")) # for supplement


  
# (2) Calculate adjacency matrix, TOM, clustering ------------------------------------


doParallel::registerDoParallel()

## Select sft
soft_power <- 3


# Calculate co-expression similarity/adjacency
set.seed(20240306)
adjacency <- adjacency(m_vsd_regress, type = "signed hybrid", power = soft_power)


## Compute topological overlap matrix (TOM)
set.seed(20240306)
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
diss_TOM <- 1 - TOM


## Run hierarchical clustering on TOM  
gene_tree <- hclust(as.dist(diss_TOM))

  
  
# (3) Module assignment -------------------------------------------------------


## Assign genes to modules across minimum module size & cut height parameter space
df_combos <- expand_grid(
    min_size = c(30, 35, 40), #c(30, 40, 50),
    cut_height = seq(from = 0.95, to = 0.98, by = 0.003)
)

set.seed(20240306)
df_modules <- map2_dfr(.x = df_combos$min_size,
                       .y = df_combos$cut_height,
                       .f = ~ tibble(sft = soft_power,
                                     min_size = .x,
                                     cut_height = .y,
                                     ensembl_gene_id = colnames(m_vsd_regress),
                                     module = cutreeDynamic(dendro = gene_tree,
                                                            distM = diss_TOM,
                                                            pamRespectsDendro = FALSE,
                                                            minClusterSize = .x,
                                                            cutHeight = .y)
                       )
) %>% 
    mutate(module = factor(module, levels = 0:(max(module))))

df_modules %>% dplyr::count(min_size, cut_height, module) %>% print(n = nrow(.))


## Plot module sizes
df_module_sizes <- df_modules %>%
    group_by(sft, min_size, cut_height) %>%
    dplyr::count(module)

df_module_sizes %>% filter(cut_height != 1) %>% 
    ggplot(aes(x = module, y = n)) +
    geom_col(width = 0.5) +
    geom_point(shape = 21, size = 0.5) +
    geom_vline(aes(xintercept = 21), lty = 2, color = "red") +
    scale_x_discrete(breaks = seq(0, 50, 10)) +
    facet_grid(min_size ~ cut_height) +
    ggtitle(paste0("SFT", soft_power, " module assignments"))


# SELECTION: SFT = 3, minimum module size = 40, cut_height = 0.98


# Assign module colors ----------------------------------------------------

## Assign colors to modules (using a different color palette than the WGCNA default)


## Select palette
igv_palette <- paletteer_d("ggsci::default_igv")


## Remove light colors and gray
igv_palette <- igv_palette[!(igv_palette %in% c("#F0E685FF", "#CDDEB7FF", "#5A655EFF", "#A9A9A9FF" ))] # n = 47


## Assign to modules in the tibble
df_colors <- tibble(
    module = 0:max((df_modules$module) %>% as.character %>% as.numeric) %>% as.factor()
) %>% 
    mutate(
        color = ifelse(module == 0, "#A9A9A9FF", igv_palette)
    )

df_modules <- df_modules %>% 
    left_join(df_colors, by = join_by(module)) %>% 
    arrange(sft, min_size, cut_height, module)


# Save final module assignments -------------------------------------------


## Save
save(df_modules, file = paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS"))

