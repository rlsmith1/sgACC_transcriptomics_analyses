
# !!!!!! RARE !!!!!!

########################################################################################

# Rare transcript level GRCCA results from Agoston's toolkit

########################################################################################

# libraries ---------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(janitor)
library(RColorBrewer)
library(biomaRt)
library(tidymodels)  
library(ggrepel)
library(raveio)
library(ggpubr)
#library(ggchicklet)


# set theme for plots -----------------------------------------------------

theme_set(theme_bw() +
            theme(plot.title = element_text(size = 11),
                  axis.title = element_text(size = 10),
                  axis.text = element_text(size = 10),
                  strip.text = element_text(size = 10),
                  legend.title = element_text(size = 8),
                  legend.text = element_text(size = 8)
            )
)

# DX COLORS
dx_colors <- c("#0072B2", "#E69F00", "#009E73", "#9966FF")
names(dx_colors) <- c("Control", "BD", "MDD", "SCZ")


# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "15Nov2023_RARE_regress_qSV1.7_ageDeath_RNAbatch_GCperc"

soft_power <- 2
minimum_size <- 35
tree_cut_height <- 0.9936

# LOAD OBJECTS 
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

# IDENTIFY MODULE SET OF INTEREST  
df_modules_filt <- df_modules %>% 
  filter(min_size == minimum_size & cut_height == tree_cut_height) %>% 
  unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
  arrange(mod_set, module)

# COVARIATES
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_clean

# DRUG MCA RESULTS
load(paste0(base_dir, "objects/drug_MCA_results.Rdata"))

# TRANSCRIPT SYMBOL TO ENSEMBL TRANSCRIPT ID MAPPING
load(paste0(base_dir, "objects/df_hsapiens_genome.RDS")) # df_hsapiens_genome

df_ensembl_to_symbol <- df_hsapiens_genome %>% 
  dplyr::select(transcript_id, transcript_name) %>% 
  distinct() %>% 
  dplyr::rename("transcript_symbol" = "transcript_name", "ensembl_transcript_id" = "transcript_id")

df_gene_to_transcript <- df_hsapiens_genome %>% 
  dplyr::select(transcript_id, transcript_name, gene_id, gene_name) %>% 
  distinct() %>% 
  dplyr::rename("gene_symbol" = "gene_name", "transcript_symbol" = "transcript_name", 
                "ensembl_gene_id" = "gene_id", "ensembl_transcript_id" = "transcript_id"
  ) %>% 
  filter(!is.na(ensembl_transcript_id))

df_modules_filt <- df_modules_filt %>% 
  left_join(df_gene_to_transcript) %>% 
  dplyr::select(ensembl_transcript_id, ensembl_gene_id, gene_symbol, module, color) %>% 
  distinct()




# define matrices ---------------------------------------------------------------

  df_order <- df_modules_filt %>% 
    dplyr::select(ensembl_transcript_id, module) %>% 
    distinct() %>% 
    arrange(module, ensembl_transcript_id)  

# X grouping vector    
  x_group <- df_order %>% 
    pull(module) %>% 
    as.numeric - 1

# X 
  X.mat <- df_vsd_regress %>% 
    dplyr::select(sample, all_of(df_order$ensembl_transcript_id)) %>% 
    column_to_rownames("sample")
  
  sample_order <- rownames(X.mat)
  
# Y and C (confound) matrices
  n_mca_dim <- 9
  YC.mat <- df_covariates_numeric %>% 
    
    # drugs to keep
    # dplyr::select(sample, age_death, brain_weight, opioids,
    #               sedative_hypnotic_anxiolitics, antidepressants, antipsychotics, mood_stabilizers, benzos) %>%
    # mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
    
    # Add MCA loadings
    left_join(df_ind_loadings, by = join_by(sample)) %>%
    dplyr::select(sample, age_death, brain_weight, suicide, all_of(paste0("dim", seq(1, n_mca_dim)))) %>%
    
    # create dx columns
    mutate(control = ifelse(grepl("control", sample), 1, 0),
           bipolar = ifelse(grepl("bipolar", sample), 1, 0),
           mdd = ifelse(grepl("mdd", sample), 1, 0),
           schizo = ifelse(grepl("schizo", sample), 1, 0)
    ) %>% 
    dplyr::select(sample, bipolar, mdd, schizo, everything(), -control) %>% 
    mutate(sample = factor(sample, levels = sample_order)) %>% 
    arrange(sample) %>% 
    column_to_rownames("sample")
  
  Y.mat <- YC.mat %>% dplyr::select(-c(brain_weight, age_death))
  head(Y.mat)
  C.mat <- YC.mat %>% dplyr::select(brain_weight, age_death)
  head(C.mat)

# CHECK ALIGNMENT
  all(rownames(X.mat) == rownames(YC.mat))
  all(colnames(X.mat) == c(df_order %>% pull(ensembl_transcript_id)))

# EXPORT FOR AGOSTON'S TOOLBOX
  type <- "RARE/"
  project_dir <- paste0(n_mca_dim, "MCA_suicide_regressAgeDeathBrainWeight/")
  cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)
  dir.create(cca_dir, recursive = TRUE)
  cca_data_dir <- paste0(cca_dir, "data/")
  dir.create(cca_data_dir)
  
  write.table(X.mat, file = paste0(cca_data_dir, "X.txt"))
  write.table(Y.mat, file = paste0(cca_data_dir, "Y.txt"))
  write.table(C.mat, file = paste0(cca_data_dir, "C.txt"))
  write.table(x_group, file = paste0(cca_data_dir, "XGroup.txt"))

# GRCCA labels files
  df_labels_x <- tibble(ensembl_transcript_id = colnames(X.mat)) %>% 
    left_join(df_modules_filt) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(Label = row_number(), .before = 1) %>% 
    dplyr::select(Label, module, ensembl_transcript_id, transcript_symbol) %>% 
    dplyr::rename_all(~c("Label", "Category", "Ensembl_id", "HGNC_id")) %>% 
    mutate(HGNC_id = ifelse(is.na(HGNC_id), Ensembl_id, HGNC_id))
  
  df_labels_y <- tibble(Label = colnames(Y.mat))
  
  write.csv(df_labels_x, paste0(cca_data_dir, "LabelsX.csv"), row.names = FALSE)
  write.csv(df_labels_y, paste0(cca_data_dir, "LabelsY.csv"), row.names = FALSE)


# GRCCA results -----------------------------------------------------------

# GRCCA RES  
  n_mca_dim <- 8
  type <- "RARE/"
  project_dir <- paste0(n_mca_dim, "MCA_suicide_regressAgeDeathBrainWeight/")
  cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)
  #analysis_dir <- "grcca_holdout1-0.20_subsamp5-0.20"
  analysis_dir <- "grcca_permutation_VARx0.1_1_mu0.1_lambda_0.9999"
  model_1 <- read_mat(paste0(cca_dir, "framework/", analysis_dir, "/res/level1/model_1.mat"))
  boot_1 <- read_mat(paste0(cca_dir, "framework/", analysis_dir, "/boot/level1/allboot_1.mat"))

# X & Y MATRICES  
Y.mat <- read.table(paste0(cca_dir, "data/Y.txt"))
df_y <- Y.mat %>% 
  rownames_to_column("sample") %>% 
  as_tibble()
covariates <- df_y %>% dplyr::select(-sample) %>% colnames

X.mat <- read.table(paste0(cca_dir, "data/X.txt"))
df_x <- X.mat %>% 
  t() %>% 
  as.data.frame %>% 
  rownames_to_column("transcript_id") %>% 
  as_tibble()
transcripts <- df_x %>% pull(transcript_id)

# READ RESULTS TABLE TO SELECT SPLIT
df_results <- read_table(paste0(cca_dir, "framework/", analysis_dir, "/res/level1/results_table.txt"))
best_split <- df_results %>% arrange(pval, -correl) %>% pull(set) %>% .[[1]]

# MODEL WEIGHTS AND STANDARD DEVIATIONS

# Y
df_y_weights <- model_1$wY %>% 
  as.data.frame %>% 
  as_tibble() %>% 
  .[best_split,] %>% 
  pivot_longer(1:ncol(.), names_to = "covariate", values_to = "weight") %>% 
  mutate(covariate = covariates)

df_y_sd <- boot_1$wY[best_split, , ] %>% 
  as.data.frame %>% 
  as_tibble() %>% 
  rename_all(~ covariates) %>% 
  map_dfr( ~ sd(.x)) %>% 
  pivot_longer(1:ncol(.), names_to = "covariate", values_to = "sd")

df_y_res <- df_y_weights %>% 
  left_join(df_y_sd) %>% 
  mutate(z_score = weight/sd) %>% 
  mutate(covariate = case_when(
    covariate == "control" ~ "Control",
    covariate == "bipolar" ~ "BD",
    covariate == "mdd" ~ "MDD",
    covariate == "schizo" ~ "SCZ",
    str_detect(covariate, "dim") ~ str_replace(covariate, "dim", "Dim "),
    TRUE ~ covariate
  )
  ) %>% 
  arrange(-abs(z_score)) %>% 
  mutate(significant = ifelse(abs(z_score) > 2, 1, 0))

# X
df_x_weights <- model_1$wX %>% 
  as.data.frame %>% 
  as_tibble() %>% 
  .[best_split,] %>% 
  pivot_longer(1:ncol(.), names_to = "ensembl_transcript_id", values_to = "weight") %>% 
  mutate(ensembl_transcript_id = transcripts)

df_x_sd <- boot_1$wX[best_split, , ] %>% 
  as.data.frame %>% 
  as_tibble() %>% 
  rename_all(~ transcripts) %>% 
  map_dfr( ~ sd(.x)) %>% 
  pivot_longer(1:ncol(.), names_to = "ensembl_transcript_id", values_to = "sd")

df_x_res <- df_x_weights %>% 
  left_join(df_x_sd) %>% 
  mutate(z_score = weight/sd) %>% 
  left_join(df_ensembl_to_symbol) %>% 
  mutate(transcript_symbol = ifelse(is.na(transcript_symbol), ensembl_transcript_id, transcript_symbol)) %>% 
  left_join(df_modules_filt) %>% 
  dplyr::select(module, color, ensembl_transcript_id, transcript_symbol, weight, sd, z_score) %>% 
  arrange(module, -abs(z_score)) %>% 
  mutate(significant = ifelse(abs(z_score) > 2, 1, 0))


# y weights ---------------------------------------------------------------

df_y_res %>% 
  mutate(significant = ifelse(abs(z_score) > 2, "*", "")) %>% 
  
  ggplot(aes(x = reorder(covariate, weight), y = weight)) +
  geom_col(aes(fill = covariate, alpha = abs(weight)), color = "black") +
  geom_errorbar(aes(ymin = weight - sd, ymax = weight + sd), width = 0.2) +
  geom_text(aes(label = significant), size = 10, vjust = 1) +
  geom_hline(aes(yintercept = 0), color = "black") +
  scale_fill_manual(values = c(dx_colors, rep("gray", 6))) +
  scale_alpha_continuous(range = c(0.1, 1)) +
  scale_x_discrete(labels = function(x) {str_replace(x, "_", " ")}) +
  ylim(c(-0.15, 0.15)) +
  labs(x = "", y = "weight", title = "Y coefficient weights",
       caption = "* = abs(z) > 2; \nError bar shows +/- 1 sd") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 0.9))

# SAVE
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/genes/genes_yres", .x),
                width = 3.5, height = 3)
)


# x weights (all) ---------------------------------------------------------------

# PLOT X RES
module_lines <- df_modules_filt %>% 
  arrange(module) %>% 
  dplyr::count(module) %>% 
  mutate(cumsum = cumsum(n)) %>% 
  pull(cumsum) + 0.5

p_x_res <- df_x_res %>% 
  arrange(module) %>% 
  mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = .$ensembl_gene_id)) %>% 
  
  ggplot(aes(x = ensembl_gene_id, y = weight)) +
  geom_col(aes(fill = I(color))) +
  geom_vline(xintercept = module_lines, color = "black", linewidth = 0.25) +
  labs(title = "X coefficient weights", x = "Gene") +
  guides(fill = guide_legend(nrow = 2)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

p_x_res

# SAVE for supplement
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/genes/x_coefficient_weights_all", .x),
                width = 10, height = 3)
)

## SIGNIFICANT ONLY
p_x_res_sig <- df_x_res %>% 
  mutate(weight = ifelse(significant == 1, weight, 0)) %>% 
  arrange(module) %>% 
  mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = .$ensembl_gene_id)) %>% 
  
  ggplot(aes(x = ensembl_gene_id, y = weight)) +
  geom_col(aes(fill = I(color))) +
  geom_vline(xintercept = module_lines, color = "black", linewidth = 0.25) +
  labs(title = "Significant X coefficient weights", x = "Gene") +
  guides(fill = guide_legend(nrow = 2)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

p_x_res_sig

# SAVE for supplement
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/genes/x_coefficient_weights_sig", .x),
                width = 10, height = 3)
)


# X modules ---------------------------------------------------------------


# HYPERGEOMETRIC TEST FOR SIGNIFICANT OVERLAP
transcript_universe <- df_x_res %>% 
  pull(ensembl_transcript_id)
grcca_transcripts <- df_x_res %>% 
  filter(abs(z_score) > 2) %>% 
  pull(ensembl_transcript_id)

df_module_overlap <- tibble()
for (mod in unique(df_x_res$module))  {
  
  mod_transcripts <- df_x_res %>% 
    filter(module == mod) %>% 
    pull(ensembl_transcript_id)
  
  mod_grcca_intersect <- intersect(grcca_transcripts, mod_transcripts)
  
  p <- 1 - phyper(q = length(mod_grcca_intersect), 
                  m = length(grcca_transcripts), 
                  n = length(gene_universe) - length(grcca_transcripts), 
                  k = length(mod_transcripts)
  )
  
  df_tmp <- tibble(module = mod,
                   overlap = length(mod_grcca_intersect),
                   mod_size = length(mod_transcripts),
                   p_value = p)
  
  df_module_overlap <- df_module_overlap %>% bind_rows(df_tmp)
  
}
df_module_overlap <- df_module_overlap %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>% 
  arrange(p_value) %>% 
  left_join(df_modules_filt %>% dplyr::select(module, color) %>% distinct()) %>% 
  mutate(module = factor(module, levels = unique(df_modules_filt$module)))

# PLOT
df_module_overlap %>%
  ggplot(aes(x = module, y = -log10(p_adj))) +
  geom_col(width = 0.05, fill = "black") +
  geom_point(aes(fill = I(color), size = overlap/mod_size), shape = 21) +
  geom_hline(aes(yintercept = -log10(0.05)), lty = 2) +
  labs(x = "module", y = "-log10(hypergeometric p-value)",
       title = "Module contribution to GRCCA gene list by hypergeometric test") +
  theme(legend.position = "none")

# SAVE
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/genes/xres_phyper_modules", .x),
                width = 6, height = 4)
)


# top weighted genes ------------------------------------------------------

df_x_res %>% 
  filter(abs(z_score) > 5) %>% 
  #arrange(-abs(weight)) %>% 
  #top_n(n = 50, wt = abs(weight)) %>% 
  
  ggplot(aes(x = reorder(transcript_symbol, weight), y = weight)) +
  geom_col(aes(fill = I(color)), color = "black") +
  geom_errorbar(aes(ymin = weight - sd, ymax = weight + sd), width = 0.1) +
  labs(title = "High confidence genes \n(abs(z-score) > 3)", x = "", y = "Weight") +
  guides(fill = guide_legend(nrow = 2)) +
  coord_flip() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "none")  

map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/genes/top_genes_abs(z)3", .x),
                width = 5, height = 10)
)

