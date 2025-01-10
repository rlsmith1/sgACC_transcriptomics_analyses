
########################################################################################

# Define modules based on control data only

########################################################################################


# libraries ---------------------------------------------------------------

  library(tidyverse)
  library(WGCNA)
  library(patchwork)
  library(janitor)
  library(RColorBrewer)


# set theme for plots -----------------------------------------------------

  theme_set(theme_bw() +
              theme(plot.title = element_text(size = 18),
                    axis.title = element_text(size = 15),
                    axis.text = element_text(size = 12),
                    strip.text = element_text(size = 15),
                    legend.title = element_text(size = 15),
                    legend.text = element_text(size = 12)))
  


# data --------------------------------------------------------------------

  prefix <- "20221209_185samples_19kgenes_vsd_qSVA17_resids"
  ctrl_prefix <- "20221209_55samplesCTRL_19kgenes_vsd_qSVA17_resids"

# LOAD REGRESSED DATA
  load(paste0("objects/", prefix, ".RDS")) # df_vsd_regress

# FILTER FOR CONTROL ONLY
  df_control <- df_vsd_regress %>% 
    filter(grepl("control", sample))

# CONVERT REGRESSED VSD COUNTS TO TRANSPOSED MATRIX
  m_control_t <- df_control %>% 
    pivot_wider(id_cols = sample, 
                names_from = ensembl_gene_id, 
                values_from = resids) %>% 
    as.data.frame %>% column_to_rownames("sample")
  

# select soft thresholding power ------------------------------------------


# CHOOSE A SET OF SOFT-THRESHOLDING POWERS TO TEST
  powers <- 1:16

# CALL THE NETWORK TOPOLOGY ANALYSIS FUNCTION
  sft <- pickSoftThreshold(m_control_t, powerVector = powers, verbose = 5)
  
  # convert results to tibble
  df_sft <- sft$fitIndices %>% as_tibble %>% clean_names
  
  # save
  save(df_sft, file = paste0("objects/", ctrl_prefix, "_SFT.RDS"))

# PLOT RESULTS   

  load(paste0("objects/", ctrl_prefix, "_SFT.RDS"))
  
  # scale-free topology fit index as a function of the soft-thresholding power
  p1 <- df_sft %>% 
    ggplot(aes(x = power, y = -sign(slope)*sft_r_sq)) +
    geom_text(aes(label = power), size = 5) +
    
    # geom_hline(yintercept = 0.8, lty = 2, color = "red") +
    geom_hline(yintercept = 0.9, lty = 2, color = "red") +
    
    xlab("Soft threshold (power)") +
    ylab("Scale free topology model fit (signed R^2)") 
  
  # Median connectivity as a function of the soft-thresholding power
  p2 <- df_sft %>% 
    ggplot(aes(x = power, y = median_k)) +
    geom_text(aes(label = power), size = 5) +
    geom_hline(yintercept = 100, lty = 2, color = "black") +
    
    xlab("Soft threshold (power)") +
    ylab("Median connectivity") # check that mean connectivity remains reasonably high (in the 100s) or above
  
  p1 + p2
  


# calculate adjacency matrix, TOM, clustering ------------------------------------

  doParallel::registerDoParallel()

# FORMAT MATRIX TO GET NAMES
  rownames(m_control_t) <- df_control %>% pull(sample) %>% unique
  colnames(m_control_t) <- df_control %>% pull(ensembl_gene_id) %>% unique

# SELECT SFT
  soft_power <- 3

# CALCULATE CO-EXPRESSION SIMILARITY AND ADJACENCY
  adjacency <- adjacency(m_control_t, power = soft_power)

# TOPOLOGICAL OVERLAP MATRIX (TOM)
  TOM <- TOMsimilarity(adjacency)
  diss_TOM <- 1 - TOM
  
  save(diss_TOM, file = paste0("objects/", ctrl_prefix, "_SFT3dissTOM.RDS"))

# HIERARCHICAL CLUSTERING ON TOM    
  gene_tree <- hclust(as.dist(diss_TOM))
  
  # plot the resulting clustering tree
  par(mfrow = c(1, 1))
  plot(gene_tree, labels = FALSE, xlab = "")



# module assignment -------------------------------------------------------


  doParallel::registerDoParallel()

# ASSIGN GENES TO MODULES AT DIFFERENT MINIMUM MODULE SIZE AND CUT HEIGHTS

  df_combos <- expand_grid(
    min_size = c(30, 50, 75),
    cut_height = seq(from = 0.95, to = 0.999, by = 0.01)
  )
  
  df_modules <- map2_dfr(.x = df_combos$min_size,
                         .y = df_combos$cut_height,
                         .f = ~ tibble(sft = 3,
                                       min_size = .x,
                                       cut_height = .y,
                                       ensembl_gene_id = df_control %>% pull(ensembl_gene_id) %>% unique,
                                       module = cutreeDynamic(dendro = gene_tree,
                                                              distM = diss_TOM,
                                                              pamRespectsDendro = FALSE,
                                                              minClusterSize = .x,
                                                              cutHeight = .y)
                         )
  )
  
  
  df_modules %>% dplyr::count(min_size, cut_height, module) %>% print(n = nrow(.))

# SAVE
  save(df_modules, file = paste0("objects/", ctrl_prefix, "_SFT3_MODS.RDS"))

  
# PLOT MODULE SIZES
  
  load(paste0("objects/", ctrl_prefix, "_SFT3_MODS.RDS")) # df_modules

  df_module_sizes <- df_modules %>%
    group_by(sft, min_size, cut_height) %>%
    dplyr::count(module)
  
  df_module_sizes %>%
    ggplot(aes(x = module, y = n)) +
    geom_point(shape = 1, size = 2) +
    facet_grid(min_size ~ cut_height) +
    ggtitle("CONTROL ONLY SFT3 module assignments")



