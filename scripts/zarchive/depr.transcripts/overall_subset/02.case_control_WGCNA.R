
##### !!!!! SUBSET !!!!

##############################################################################################################

# OVERALL TRANSCRIPTS: Select soft-thresholding power, calculate adjacency & TOM, cluster, and assign modules

##############################################################################################################

# libraries ---------------------------------------------------------------

library(tidyverse)
library(WGCNA)
library(patchwork)
library(janitor)
library(RColorBrewer)
library(ggpubr)
library(paletteer)
library(segmented)



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

# LOAD GENE NAMES
load(paste0(base_dir, "objects/08Mar2024_GENES_qSVAgeSexRaceGC.RDS")) # df_vsd_regress (gene_level)
gene_ids <- df_vsd_regress %>% dplyr::select(-sample) %>% colnames
length(gene_ids) # n = 18,677
rm(list = "df_vsd_regress")

# LOAD TRANSCRIPT DATA
prefix <- "08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC"

# RAW TRANSCRIPT COUNTS
df_transcript_raw_counts <- read.csv(paste0(base_dir, "data/transcript_count_matrix.csv")) %>% 
  as_tibble() %>% 
  dplyr::rename("ensembl_transcript_id" = "X") %>% 
  rename_all(~str_remove(.x, "_MergedBam_stringtieOutput")) %>% 
  clean_names() %>% 
  rename_all(~str_remove(.x, "x")) 

# LOAD REGRESSED DATA
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress

# HSAPIENS GENOME
load(paste0(base_dir, "objects/hsapiens_genome_v110.RDS")) # df_hsapiens_genome


# OPTIONAL? map transcripts to genes used in gene-level analysis ---------------------------------------------


ncol(df_vsd_regress) - 1 # n = 69769 (low count transcripts were removed in 01.raw_count_preprocessing.R)

## FILTER FOR TRANSCRIPTS THAT ARE ASSOCIATED WITH GENES WE INCLUDED IN ANALYSIS (SUBS)

# map genes to transcripts
df_gene_to_transcript <- df_hsapiens_genome %>% 
  dplyr::select(gene_id, transcript_id) %>% 
  distinct() %>% 
  filter(!is.na(transcript_id))

# filter for genes that were involved in the gene-level analysis
transcript_ids <- df_gene_to_transcript %>% 
  filter(gene_id %in% gene_ids) %>% 
  pull(transcript_id) %>% 
  unique
length(transcript_ids) # map to 158115 transcripts

# select columns of transcripts involved in this analysis
df_vsd_regress_all <- df_vsd_regress %>% 
  dplyr::select(sample, any_of(transcript_ids)) 
ncol(df_vsd_regress_all) - 1 # n = 62,773

# load(paste0(base_dir, "objects/22Nov2023_OVERALLsubs_regress_qSV1.7_ageDeath_RNAbatch_GCperc.RDS")) # df_vsd_regress_all # n = 52327
# load(paste0(base_dir, "objects/08Feb2024_OVERALLsubs_regress_qSV1.7_ageDeath_RNAbatch_GCperc.RDS")) # df_vsd_regress_all # n = 62371
# load(paste0(base_dir, "objects/10Feb2024_OVERALL_regress_qSV1.7_ageDeath_RNAbatch_GCperc.RDS")) # df_vsd_regress_all # n = 69116

# # SAVE  
# prefix2 <- "05Mar2024_OVERALLsubs_regress_qSV1.7_ageDeath_RNAbatch_GCperc"
# save(df_vsd_regress_all, file = paste0(base_dir, "objects/", prefix2, ".RDS"))
# 
# # SAVE AS A MATRIX FOR THE CLUSTER
# count_matrix <- df_vsd_regress_all %>% column_to_rownames("sample")
# save(count_matrix, file = paste0(base_dir, "objects/", prefix2, "_matrix.RDS"))


# OPTIONAL? Additional filtering: by coefficient of variance ------------------------

## "Selecting information features and reducing dimensionality"
# https://www.nature.com/articles/s41576-023-00586-w
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-193 (pick arbitrary threshold?)

# However, we don't just want to filter on variance because variance scales with mean and some transcripts have very low expression
# So we can use coef of variance, which filters on variance relative to mean


# FILTER BY MEAN/SD (COEF OF VAR) RAW COUNTS
df_mean_sd <- df_transcript_raw_counts %>% 
  filter(ensembl_transcript_id %in% colnames(df_vsd_regress)) %>% 
  pivot_longer(2:ncol(.), names_to = "sample", values_to = "raw_counts") %>% 
  group_by(ensembl_transcript_id) %>% 
  summarise(
    mean = mean(raw_counts, na.rm = TRUE),
    stdev = sd(raw_counts, na.rm = TRUE),
    CV = stdev/mean
  ) %>% 
  arrange(-CV)
df_mean_sd <- df_mean_sd %>% mutate(z_CV = scale(CV)[,1])

# # FILTER BY MEAN/SD (COEF OF VAR) NORMALIZED/REGRESSED COUNTS
# df_mean_sd <- df_vsd_regress %>% 
#   pivot_longer(2:ncol(.), names_to = "ensembl_transcript_id", values_to = "residuals") %>% 
#   group_by(ensembl_transcript_id) %>% 
#   summarise(
#     mean = mean(residuals, na.rm = TRUE),
#     stdev = sd(residuals, na.rm = TRUE),
#     CV = stdev/mean
#   ) %>% 
#   arrange(-CV)
# df_mean_sd <- df_mean_sd %>% mutate(z_CV = scale(CV)[,1])


# PLOTS

# plot mean vs sd
df_mean_sd %>% 
  ggplot(aes(x = mean %>% log10, y = stdev %>% log10)) +
  geom_point(shape = 21) +
  geom_density2d() +
  geom_smooth(method = "lm", color = "red", lty = 2) +
  stat_cor(label.sep = "\n") +
  labs(x = "log10(mean)", y = "log10(standard deviation)",
       title = "Relationship between mean and standard deviation in transcript raw counts")
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/mean_variance_relationship", .x), width = 6, height = 4)
)

# plot mean vs CV
df_mean_sd %>% 
  ggplot(aes(x = log10(mean), y = CV %>% log10)) +
  geom_point(shape = 21) +
  geom_density2d() +
  labs(x = "log10(mean)", y = "log10(coefficient of variation)",
       title = "Relationship between mean and coefficient of variation in transcript raw counts")
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/mean_CV_relationship", .x), width = 6, height = 4)
)


# by z-score?
df_mean_sd %>% filter(z_CV < -1)
df_mean_sd %>% filter(CV < 0.25)
df_mean_sd %>% arrange(mean)
summary(df_mean_sd %>% pull(mean))


# plot
quantiles <- summary(df_mean_sd %>% pull(CV))

cutoff <- quantiles[2]
df_mean_sd %>% 
  ggplot(aes(x = CV)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(y = 3.0), position = position_jitter(0.001), alpha = 0.5, shape = 21) +
  geom_density(color = "blue", fill = "transparent") +
  geom_boxplot(aes(y = -0.1), color = "blue", width = 0.2, outlier.shape = NA) +
  geom_vline(xintercept = cutoff, lty = 2, color = "red") +
  labs(x = "coefficient of variation = standard deviation / mean", y = "",
       title = "Coefficient of variation across transcripts",
       caption = "Red dashed line indicates coefficient of variation cutoff (0.36)") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  ) 
# low CV means the transcript expression doesn't really change across samples, thus it will not be informative in our GRCCA analysis (or WGCNA for that matter)

map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/coefficient_of_variation", .x), width = 6, height = 4)
)

# keep transcripts that are above CV cutoff
keep <- df_mean_sd %>% filter(CV >= cutoff) %>% pull(ensembl_transcript_id)
length(keep) # n = 54302 (CVq1; where q1 ~= 0.36)

df_vsd_regress_filt <- df_vsd_regress %>% dplyr::select(sample, any_of(keep)) 
ncol(df_vsd_regress_filt) - 1 # n = 54302

# SAVE  
prefix2 <- paste0(prefix, "_CVq1")
save(df_vsd_regress_filt, file = paste0(base_dir, "objects/", prefix2, ".RDS")) # df_vsd_regress_filt

# SAVE AS A MATRIX FOR THE CLUSTER
count_matrix <- df_vsd_regress_filt %>% column_to_rownames("sample")
save(count_matrix, file = paste0(base_dir, "objects/", prefix2, "_matrix.RDS"))



# select soft thresholding power ------------------------------------------

####### ****************** RUN ON CLUSTER *********************************

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


# assign module colors & resave ----------------------------------------------------

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
  
  write.csv(paste0(base_dir, "outputs/tables/for_manuscript/", prefix2, "_SIGNED_SFT", soft_power, "_MODS.csv"),
            row.names = FALSE)


