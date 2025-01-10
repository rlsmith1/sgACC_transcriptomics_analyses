
######################################################################################

### TRANSCRIPTS
## Read in raw count data, normalize using DESeq2, qSVA correction & regress covariates

######################################################################################

# libraries ---------------------------------------------------------------

library(tidyverse)
library(janitor)
library(RColorBrewer)
library(DESeq2)
library(patchwork)
library(biomaRt)
library(corrr)
library(ggrepel)


# Load transcript data ---------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"

# COVARIATES
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)

# RAW TRANSCRIPT COUNTS
df_transcript_raw_counts <- read.csv(paste0(base_dir, "data/transcript_count_matrix.csv")) %>% 
  as_tibble() %>% 
  dplyr::rename("ensembl_transcript_id" = "X") %>% 
  rename_all(~str_remove(.x, "_MergedBam_stringtieOutput")) %>% 
  clean_names() %>% 
  rename_all(~str_remove(.x, "x")) 

# filter out samples that weren't in genes
df_transcript_raw_counts <- df_transcript_raw_counts %>% dplyr::select(ensembl_transcript_id, contains(df_covariates$sample))

# LOAD QSVS
load(paste0(base_dir, "objects/06March2024_qSVs.Rdata"))


# COMBINE qSVS WITH COVARIATE DATA  
df_covar_qsv <- df_covariates_numeric %>% left_join(df_qsvs)


# Remove low count transcripts --------------------------------------------

# At least 10 raw counts in at least 80% of samples
keep <- df_transcript_raw_counts %>% 
  mutate_if(is.numeric, ~ifelse(.x >= 10, 1, 0)) %>% 
  mutate(thresh = rowSums(dplyr::select(., -ensembl_transcript_id))/(ncol(.) - 1), .before = 2) %>% 
  filter(thresh >= 0.80) %>% 
  pull(ensembl_transcript_id)
length(keep) # 72403

# filter df_transcripts for these
df_transcript_raw_counts_filtered <- df_transcript_raw_counts %>% 
  filter(ensembl_transcript_id %in% keep)


# normalize transcript counts using VST from DESEq2 ---------------------------------------------

m_counts <- df_transcript_raw_counts_filtered %>% 
  as.data.frame %>% 
  column_to_rownames("ensembl_transcript_id") %>% 
  as.matrix + 1 # add a pseudo-count of 1 to avoid zeros in count mat

m_vsd <- varianceStabilizingTransformation(m_counts, blind = TRUE, fitType = "parametric")


# PCA1: on normalized counts with no regression ---------------------------------------------------------------------


# PCA ON NORMALIZED COUNTS 
pca1 <- prcomp(m_vsd,
               center = TRUE, scale = TRUE)

df_pca <- pca1$rotation %>% 
  as.data.frame %>% 
  rownames_to_column("sample") %>% 
  as_tibble() %>% 
  left_join(df_covar_qsv) %>% 
  clean_names()

df_var <- summary(pca1) %>% 
  .$importance %>% 
  as.data.frame %>% 
  rownames_to_column("metric") %>% 
  as_tibble() %>% 
  pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
  mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
  pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
  clean_names()

# PLOT PCA
p_pca <- df_pca %>% 
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(fill = q_sv1), shape = 21, size = 2) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  guides(fill = guide_colorbar(title = "qSV1")) +
  
  # add origin lines
  # geom_vline(xintercept = 0, color = "black", lty = 2) +
  # geom_hline(yintercept = 0, color = "black", lty = 2) +
  
  # give % variance in axis labels
  xlab(paste0("PC1: ", 
              df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
              "% of variance")) +
  ylab(paste0("PC2: ", 
              df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
              "% of variance")) +
  theme(legend.position = c(0.10, 0.80))

# plot eigenvalues
p_var <- df_var %>% 
  ggplot(aes(x = pc, y = proportion_of_variance)) +
  geom_bar(aes(fill = proportion_of_variance), stat = "identity",
           color = "black", linewidth = 0.05) +
  geom_hline(aes(yintercept = 0.02), lty = 2, color = "red") +
  scale_fill_gradient(low = "white", high = "midnightblue", 
                      limits = c(0, 1.0), guide = "none") +
  labs(x = "PC", y = "Proportion of variance explained")

p_pca + p_var + plot_annotation(title = "PCA 1: VST normalized counts before qSV regression")

# SAVE
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/qSVA/PCA1_TRANSCRIPTS_normalized_counts", .x),
                width = 10, height = 5)
)

# CORRELATE PCs WITH COVARIATES
df_pc_covar_cor <- df_pca %>% 
  pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
  pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
  group_by(pc, covariate) %>% 
  nest() %>% 
  mutate(pearsons_r = map(.x = data, .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate),
         p_val = map(.x = data, .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value)) %>% 
  unnest(cols = c(pearsons_r, p_val)) %>% 
  ungroup %>% 
  mutate(p_adj = p.adjust(p_val, method = "fdr"))

# PLOT COVARIATE PC CORRELATION

# set order
m_covariates <- df_covariates_numeric %>% column_to_rownames("sample") %>% t()
covariate_order <- m_covariates[hclust(dist(m_covariates))$order,] %>% rownames

# plot
df_pc_covar_cor %>% 
  mutate(covariate = str_replace(covariate, "q_sv", "qSV"),
         covariate = factor(covariate, levels = c(covariate_order, paste0("qSV", 1:16))),
         pc = str_remove(pc, "pc") %>% as.numeric) %>% 
  left_join(df_var) %>% 
  filter(proportion_of_variance > 0.02) %>% 
  mutate(significant = ifelse(p_adj < 0.05, "yes", "no")) %>% 
  
  ggplot(aes(x = pc, y = covariate)) +
  geom_tile(aes(fill = pearsons_r, color = significant), width = 0.98, height = 0.98, linewidth = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  scale_color_manual(values = c("yes" = "black", "no" = "transparent"), guide = "none") +
  guides(fill = guide_colorbar(title = "Pearson's r")) +
  labs(x = "PC (> 2% explained var)", title = "PC-covariate correlations",
       caption = "Box outline indicates FDR < 0.05") +
  theme_classic() +
  theme(legend.key.height = unit(2, "cm"))

# SAVE
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/qSVA/PCA1_TRANSCRIPTS_covariate_correlations", .x),
                width = 5, height = 8)
)



# REGRESS qSVs from transcript level --------------------------------------------


# REGRESS qSVs THAT MEET CRITERIA AND AGE_DEATH (SAME AS GENE LEVEL)
df_qsvs %>% dplyr::select(-sample) %>% colnames
keep_qsvs <- paste0("qSV", c(seq(1, 7), 9)) # correlate with known covariates or var explained > 2%
mod <- df_covar_qsv %>% 
  dplyr::select(all_of(keep_qsvs), sex_at_birth, race, age_death, gc_percent) %>% 
  as.data.frame %>% 
  as.matrix

# RUN REGRESSION  
doParallel::registerDoParallel()
df_vsd_regress <- tibble(sample = colnames(m_vsd))  
for (i in 1:nrow(m_vsd)) {
  
  if (i %% 100 == 0) {print(paste0(i, " of ", nrow(m_vsd)))}
  
  transcript <- rownames(m_vsd)[i]
  df <- cbind(m_vsd = m_vsd[i,], mod) %>% as.data.frame()
  m_vsd_clean <- (lm(m_vsd ~ ., data = df) %>% residuals) #+ mean(df$m_vsd) # add mean for interpretable expression value
  
  df_vsd_regress[[transcript]] <- m_vsd_clean
  
}  

df_vsd_regress
  


# PCA2: vsd ~ qSV1-7 + 9 + age_death + sex_at_birth + race + gc_percent --------------------------------------
  
# PCA ON NORMALIZED & REGRESSED COUNTS    
  pca <- prcomp(df_vsd_regress %>% 
                  column_to_rownames("sample") %>% 
                  t(),
                center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    left_join(df_covar_qsv) %>% 
    clean_names()
  
  df_var <- summary(pca) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()
  
# PLOT PCA
  p_pca <- df_pca %>% 
    ggplot(aes(x = pc1, y = pc2)) +
    geom_point(aes(fill = q_sv1), shape = 21, size = 2) +
    scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) +
    guides(fill = guide_colorbar(title = "qSV1")) +
    
    # add origin lines
    # geom_vline(xintercept = 0, color = "black", lty = 2) +
    # geom_hline(yintercept = 0, color = "black", lty = 2) +
    
    # give % variance in axis labels
    xlab(paste0("PC1: ", 
                df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance")) +
    ylab(paste0("PC2: ", 
                df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance")) +
    theme(legend.position = c(0.90, 0.26))
  
  # plot eigenvalues
  p_var <- df_var %>% 
    ggplot(aes(x = pc, y = proportion_of_variance)) +
    geom_bar(aes(fill = proportion_of_variance), stat = "identity",
             color = "black", linewidth = 0.05) +
    geom_hline(aes(yintercept = 0.01), lty = 2, color = "red") +
    scale_fill_gradient(low = "white", high = "midnightblue", 
                        limits = c(0, 0.015), guide = "none") +
    labs(x = "PC", y = "Proportion of variance explained")
  
  p_pca / p_var + plot_annotation(title = "PCA 3: VST normalized counts with qSV + race + sex at birth + age at death + GC percent regression")
  
  # save
  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/qSVA/PCA3_TRANSCRIPTS_normalized_and_regressed_counts", .x),
                  width = 7, height = 8)
  )
  
  
# CORRELATE PCs WITH COVARIATES
  df_pc_covar_cor <- df_pca %>% 
    pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
    pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
    group_by(pc, covariate) %>% 
    nest() %>% 
    mutate(pearsons_r = map(.x = data, .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate),
           p_val = map(.x = data, .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value)) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    ungroup %>% 
    mutate(p_adj = p.adjust(p_val, method = "fdr"))
  
  # plot  
  df_pc_covar_cor %>% 
    mutate(covariate = str_replace(covariate, "q_sv", "qSV"),
           covariate = factor(covariate, levels = c(covariate_order, paste0("qSV", 1:16))),
           pc = str_remove(pc, "pc") %>% as.numeric) %>% 
    left_join(df_var) %>% 
    filter(proportion_of_variance > 0.01) %>% 
    mutate(significant = ifelse(p_adj < 0.05, "yes", "no")) %>% 
    
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = pearsons_r, color = significant), width = 0.98, height = 0.98, linewidth = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    scale_color_manual(values = c("yes" = "black", "no" = "transparent"), guide = "none") +
    guides(fill = guide_colorbar(title = "Pearson's r")) +
    labs(x = "PC (> 1% explained var)", title = "PC-covariate correlations following qSV + race + sex at birth + age at death + GC percent regression",
         caption = "Box outline indicates FDR < 0.05") +
    theme_classic() +
    theme(legend.key.height = unit(2, "cm"))
  
  # SAVE
  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/qSVA/PCA3_TRANSCRIPTS_covariate_correlations", .x),
                  width = 8, height = 8)
  )
  
  

# save object for WTCNA ---------------------------------------------------


  # SAVE OBJECT  
  prefix <- "08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC"
  save(df_vsd_regress, file = paste0(base_dir, "objects/", prefix, ".RDS"))
  
  
  
  
  