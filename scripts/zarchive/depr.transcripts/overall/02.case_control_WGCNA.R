
##### !!!!! ALL !!!!

##############################################################################################################

# OVERALL TRANSCRIPTS: Select soft-thresholding power, calculate adjacency & TOM, cluster, and assign modules

##############################################################################################################

# libraries ---------------------------------------------------------------

library(tidyverse)
library(WGCNA)
library(sva)
library(patchwork)
library(janitor)
library(RColorBrewer)


# set theme for plots -----------------------------------------------------

theme_set(theme_bw() +
            theme(plot.title = element_text(size = 12),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 10),
                  strip.text = element_text(size = 10),
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 10)))

# data --------------------------------------------------------------------



base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "15Nov2023_TRANSCRIPTS_regress_qSV1.7_ageDeath_RNAbatch_GCperc"

# RAW TRANSCRIPT COUNTS
df_transcript_raw_counts <- read.csv(paste0(base_dir, "data/transcript_count_matrix.csv")) %>% 
  as_tibble() %>% 
  dplyr::rename("ensembl_transcript_id" = "X") %>% 
  rename_all(~str_remove(.x, "_MergedBam_stringtieOutput")) %>% 
  clean_names() %>% 
  rename_all(~str_remove(.x, "x")) 

# LOAD REGRESSED DATA
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress



# filter for transcripts that have greater than 10 counts in 80% of samples ---------------------------------------------


# Greater than 10 counts but less than 100 counts in 80% of samples

keep <- df_transcript_raw_counts %>% 
  mutate_if(is.numeric, ~ifelse(.x > 10, 1, 0)) %>% 
  mutate(thresh = rowSums(dplyr::select(., -ensembl_transcript_id))/(ncol(.) - 1), .before = 2) %>% 
  filter(thresh > 0.80) %>% 
  pull(ensembl_transcript_id)

length(keep)

df_vsd_regress_all <- df_vsd_regress %>% dplyr::select(sample, any_of(keep)) # n = 69,116

# SAVE  
prefix2 <- "15Nov2023_OVERALL_regress_qSV1.7_ageDeath_RNAbatch_GCperc"
save(df_vsd_regress_all, file = paste0(base_dir, "objects/", prefix2, ".RDS")) # df_vsd_regress_all


# SAVE AS A MATRIX FOR THE CLUSTER
count_matrix <- df_vsd_regress_all %>% column_to_rownames("sample")
save(count_matrix, file = paste0(base_dir, "objects/", prefix2, "_matrix.RDS"))


# select soft thresholding power ------------------------------------------

#### CODE FOR CLUSTER

# IN LOCAL SPACE: transfer matrix to cluster
# cd ~/Documents/PhD/projects/sgacc_wgcna_grcca/objects
# rsync 10Feb2024_OVERALL_regress_qSV1.7_ageDeath_RNAbatch_GCperc_matrix.RDS smithral@biowulf.nih.gov:/data/smithral/sgacc_wgcna/transcripts/data

# ON CLUSTER
# ssh smithral@biowulf.nih.gov
# cd /data/smithral/sgacc_wgcna/transcripts
# spersist --mem=500g --gres=lscratch:800
# module load rstudio R
# Rscript scripts/00.calc_sft.R

# IN LOCAL SPACE
# cd ~/Documents/PhD/projects/sgacc_wgcna_grcca/objects
# rsync smithral@biowulf.nih.gov:/data/smithral/sgacc_wgcna/transcripts/outputs/08Feb2024_OVERALLsubs_regress_qSV1.7_ageDeath_RNAbatch_GCperc_SFT.RDS .

###

# PLOT RESULTS   

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
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/transcripts/OVERALL/sft", .x),
                width = 8, height = 4)
)


# (CLUSTER) calculate adjacency matrix, TOM, clustering ------------------------------------


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




# (CLUSTER) module assignment -------------------------------------------------------


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


# #### plot results from cluster ------------------------------------------


# PLOT MODULE SIZES
soft_power <- 2
prefix2 <- "15Nov2023_OVERALL_regress_qSV1.7_ageDeath_RNAbatch_GCperc"
load(paste0(base_dir, "objects/", prefix2, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

df_module_sizes <- df_modules %>%
  group_by(sft, min_size, cut_height) %>%
  dplyr::count(module)

df_module_sizes %>% filter(cut_height > 0.98 & cut_height < 0.99) %>% 
  ggplot(aes(x = module, y = n)) +
  geom_col(width = 0.2) +
  geom_point() +
  geom_vline(aes(xintercept = 20), lty = 2, color = "red") +
  scale_x_discrete(breaks = seq(0, 250, 10)) +
  facet_grid(min_size ~ cut_height) +
  ggtitle(paste0("OVERALL SFT", soft_power, " module sizes"))


map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/transcripts/OVERALL/sft", soft_power, "_module_sizes", .x), width = 10, height = 5)
)


