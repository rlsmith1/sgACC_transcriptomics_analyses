
##### !!!!! NOTE: transcripts have already been filtered pre-WTCNA! !!!!

#########################################################################################################################

# TRANSCRIPTS: Select soft-thresholding power, calculate adjacency & TOM, cluster, and assign modules

#########################################################################################################################

####### ****************** RUN ON CLUSTER ********************************* #######

### Setup ###
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/load_transcripts.R"))

# Set directories to save objects
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")


# (CLUSTER) Select soft thresholding power ------------------------------------------


# # CHOOSE A SET OF SOFT-THRESHOLDING POWERS TO TEST
# powers <- 1:16
# 
# # CALL THE NETWORK TOPOLOGY ANALYSIS FUNCTION
# sft <- pickSoftThreshold(m_vsd_regress_t, powerVector = powers, verbose = 5)
# 
# # convert results to tibble
# df_sft <- sft$fitIndices %>% as_tibble %>% clean_names
# 
# # save
# save(df_sft, file = paste0(base_dir, "objects/", prefix2, "_SFT.RDS"))

# PLOT RESULTS   

#### CODE FOR CLUSTER

# IN LOCAL SPACE: transfer matrix to cluster
# cd ~/Documents/PhD/projects/sgacc_wgcna_grcca/objects
# rsync 08Feb2024_OVERALLsubs_regress_qSV1.7_ageDeath_RNAbatch_GCperc_matrix.RDS smithral@biowulf.nih.gov:/data/smithral/sgacc_wgcna/transcripts/data

# ON CLUSTER
# ssh smithral@biowulf.nih.gov
# cd /data/smithral/sgacc_wgcna/transcripts
# spersist --mem=500g --gres=lscratch:800
# module load rstudio R
# Rscript scripts/00.calc_sft.R

# IN LOCAL SPACE
# cd ~/Documents/PhD/projects/sgacc_wgcna_grcca/objects
# rsync smithral@biowulf.nih.gov:/data/smithral/sgacc_wgcna/transcripts/outputs/08Feb2024_OVERALLsubs_regress_qSV1.7_ageDeath_RNAbatch_GCperc_SFT.RDS .

load(paste0(base_dir, "objects/", prefix2, "_SFT.RDS")) # (load from cluster)

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

# SAVE
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/transcripts/OVERALLsubs/sft", .x),
                width = 8, height = 4)
)


# (CLUSTER) Calculate adjacency matrix, TOM, clustering ------------------------------------


# LOAD DATA
base_dir <- "/Users/work/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/projects/WGCNA/acsg_wgcna/"
prefix2 <- "20230315_185samples_25kRAREtranscripts_vst_qSVA123567_MP_RNAbatch_Race_resids"
load(paste0(base_dir, "objects/", prefix2, ".RDS")) # df_vsd_regress_rare
load(paste0("objects/", prefix2, "_matrix.RDS")) # m_vsd_regress_t

### TOO BIG - RUN ON BIOWULF  

#Sys.setenv('R_MAX_VSIZE'=32000000000)
doParallel::registerDoParallel()

# CONVERT REGRESSED VSD COUNTS TO TRANSPOSED MATRIX
# m_vsd_regress_t <- df_vsd_regress_rare %>%
#   pivot_wider(id_cols = sample,
#               names_from = ensembl_transcript_id,
#               values_from = resids) %>%
#   as.data.frame %>% column_to_rownames("sample")
# 
# rownames(m_vsd_regress_t) <- df_vsd_regress_rare %>% pull(sample) %>% unique
# colnames(m_vsd_regress_t) <- df_vsd_regress_rare %>% pull(ensembl_transcript_id) %>% unique

# save matrix for export to Biowulf
# save(m_vsd_regress_t, file = paste0("objects/", prefix2, "_matrix.RDS")) # m_vsd_regress_t


# SELECT SFT
soft_power <- 3

# CALCULATE CO-EXPRESSION SIMILARITY AND ADJACENCY
adjacency <- adjacency(m_vsd_regress_t, type = "signed hybrid", power = soft_power)

# TOPOLOGICAL OVERLAP MATRIX (TOM)
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
rm(list = "adjacency")

diss_TOM <- 1 - TOM
rm(list = "TOM")

# HIERARCHICAL CLUSTERING ON TOM    
gene_tree <- hclust(as.dist(diss_TOM))

# plot the resulting clustering tree
par(mfrow = c(1, 1))
plot(gene_tree, labels = FALSE, xlab = "")




# (CLUSTER) Module assignment -------------------------------------------------------


# PACKAGE RECOMMENDATIONS


# ASSIGN GENES TO MODULES AT DIFFERENT MINIMUM MODULE SIZE AND CUT HEIGHTS

doParallel::registerDoParallel()

df_combos <- expand_grid(
  min_size = c(30, 40, 50),
  cut_height = seq(from = 0.998, to = 0.999, by = 0.0001)
)

df_modules <- map2_dfr(.x = df_combos$min_size,
                       .y = df_combos$cut_height,
                       .f = ~ tibble(sft = soft_power,
                                     min_size = .x,
                                     cut_height = .y,
                                     ensembl_transcript_id = unique(df_vsd_regress_rare$ensembl_transcript_id),
                                     module = cutreeDynamic(dendro = gene_tree,
                                                            distM = diss_TOM,
                                                            pamRespectsDendro = FALSE,
                                                            minClusterSize = .x,
                                                            cutHeight = .y)
                       )
) %>% 
  mutate(module = factor(module, levels = 0:(max(module))))

df_modules %>% dplyr::count(min_size, cut_height, module) %>% print(n = nrow(.))


# ASSIGN COLORS TO MODULES  

df_modules <- df_modules %>% 
  mutate(color = labels2colors(as.numeric(as.character(module)))) %>% 
  
  mutate(color = case_when(
    
    color == "yellow" ~ "orange",
    color == "lightcyan" ~ "darkcyan",
    color == "lightgreen" ~ "seagreen",
    color == "lightyellow" ~ "darkorange",
    TRUE ~ color
    
  ))

# SAVE
save(df_modules, file = paste0("objects/", prefix2, "_SIGNED_SFT", soft_power, "_MODS.RDS"))


# #### Plot results from cluster ------------------------------------------

# in local space: 
# cd /Users/smithral/Documents/PhD/projects/sgacc_wgcna_grcca/objects
# rsync smithral@biowulf.nih.gov:/data/smithral/sgacc_wgcna/transcripts/outputs/05Mar2024_OVERALLsubs_regress_qSV1.7_ageDeath_RNAbatch_GCperc_map2genes_CV_SIGNED_SFT2_MODS.RDS . 


# PLOT MODULE SIZES
soft_power <- 2
load(paste0(base_dir, "objects/", prefix2, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

# CALCULATE MODULE SIZES
df_module_sizes <- df_modules %>%
  group_by(sft, min_size, cut_height) %>%
  dplyr::count(module)

# where do we lose the gray module?
df_module_sizes %>% filter(module == 0) %>% dplyr::count(min_size, cut_height) %>% print(n = nrow(.))

# PLOT
df_module_sizes %>% #filter(cut_height != 1) %>% 
  ggplot(aes(x = module, y = n)) +
  geom_col(width = 0.5) +
  geom_point(shape = 21, size = 0.5) +
  geom_vline(aes(xintercept = 21), lty = 2, color = "red") +
  scale_x_discrete(breaks = seq(0, 50, 10)) +
  facet_grid(min_size ~ cut_height) +
  ggtitle(paste0("SFT", soft_power, " module assignments"))

# TRY SFT = 2, mininum module size = 35, cut height = 0.987

map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/sft", soft_power, "_module_sizes", .x), width = 12, height = 6)
)


# Assign module colors & resave ----------------------------------------------------

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
  dplyr::select(-color) %>% 
  left_join(df_colors, by = join_by(module)) %>% 
  arrange(sft, min_size, cut_height, module)

# SAVE
save(df_modules, file = paste0(base_dir, "objects/", prefix2, "_SIGNED_SFT", soft_power, "_MODS.RDS"))

# write to csv  
df_modules %>% 
  unite(col = "Module set", c(sft, min_size, cut_height), sep = "_") %>% 
  dplyr::select(ensembl_transcript_id, everything()) %>% 
  group_by(`Module set`) %>% 
  arrange(module) %>% 
  
  write.csv(paste0(base_dir, "outputs/tables/TableSXX", prefix2, "_SIGNED_SFT", soft_power, "_MODS.csv"),
            row.names = FALSE)


