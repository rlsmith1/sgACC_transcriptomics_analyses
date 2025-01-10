
########################################################################################

# Select soft-thresholding power, calculate adjacency & TOM, cluster, and assign modules

########################################################################################


### Set directories & load data ###

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/FigS3_WGCNAenrichment/")
tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_genes.R"))


# Select soft thresholding power ------------------------------------------

  
# CONVERT REGRESSED VSD COUNTS TO TRANSPOSED MATRIX
  m_vsd_regress <- df_vsd_regress %>% column_to_rownames("sample")
  
# CHOOSE A SET OF SOFT-THRESHOLDING POWERS TO TEST
  powers <- 1:16

# CALL THE NETWORK TOPOLOGY ANALYSIS FUNCTION
  sft <- pickSoftThreshold(m_vsd_regress, powerVector = powers, verbose = 5)
  
  # convert results to tibble
  df_sft <- sft$fitIndices %>% as_tibble %>% clean_names
  
  # save
  save(df_sft, file = paste0(base_dir, "objects/", prefix, "_SFT.RDS")) # for downstream analyses
  save(df_sft, file = paste0(analysis_objects_dir, prefix, "_SFT.RDS")) # for supplement


  
# Calculate adjacency matrix, TOM, clustering ------------------------------------

  doParallel::registerDoParallel()
  
# CONVERT REGRESSED VSD COUNTS TO TRANSPOSED MATRIX
  m_vsd_regress <- df_vsd_regress %>% column_to_rownames("sample")
  
# SELECT SFT
  soft_power <- 3

# CALCULATE CO-EXPRESSION SIMILARITY AND ADJACENCY
  set.seed(20240306)
  adjacency <- adjacency(m_vsd_regress, type = "signed hybrid", power = soft_power)

# TOPOLOGICAL OVERLAP MATRIX (TOM)
  set.seed(20240306)
  TOM <- TOMsimilarity(adjacency, TOMType = "signed")
  diss_TOM <- 1 - TOM

# HIERARCHICAL CLUSTERING ON TOM    
  gene_tree <- hclust(as.dist(diss_TOM))
  
  # plot the resulting clustering tree
  par(mfrow = c(1, 1))
  plot(gene_tree, labels = FALSE, xlab = "", main = "Gene dendrogram") #main = paste0("SFT", soft_power, " gene tree"))

# PLOT TOM
  # library(pheatmap)
  # library(viridis)
  # m_cor <- cor(m_vsd_regress_t)
  # pheatmap(m_cor[4000:5000, 4000:5000],
  #          color = viridis(10),
  #          border_color = NA,
  #          show_colnames = FALSE, 
  #          show_rownames = FALSE,
  #          treeheight_row = 0,
  #          treeheight_col = 0,
  #          main = "Gene co-expression matrix")
  
  
# Module assignment -------------------------------------------------------


# ASSIGN GENES TO MODULES AT DIFFERENT MINIMUM MODULE SIZE AND CUT HEIGHTS
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


# PLOT MODULE SIZES
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

  # map(
  #   .x = c(".png", ".pdf"),
  #   .f = ~ ggsave(paste0(figures_dir, "module_sizes", .x), width = 12, height = 5)
  # )
  
  # SELECTION: SFT = 3, minimum module size = 40, cut_height = 0.98
  

# Assign module colors ----------------------------------------------------

# ASSIGN COLORS TO MODULES 

  # select palette
  igv_palette <- paletteer_d("ggsci::default_igv")
  
  # remove light colors & gray
  igv_palette <- igv_palette[!(igv_palette %in% c("#F0E685FF", "#CDDEB7FF", "#5A655EFF", "#A9A9A9FF" ))] # n = 47
  
  # assign to modules in tibble
  df_colors <- tibble(
    module = 0:max((df_modules$module) %>% as.character %>% as.numeric) %>% as.factor()
  ) %>% 
    mutate(
      color = ifelse(module == 0, "#A9A9A9FF", igv_palette)
    )
  
  # add to df_modules tibble
  df_modules <- df_modules %>% 
    left_join(df_colors, by = join_by(module)) %>% 
    arrange(sft, min_size, cut_height, module)


# Save final module assignments -------------------------------------------

    
# SAVE
  save(df_modules, file = paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS"))
  
