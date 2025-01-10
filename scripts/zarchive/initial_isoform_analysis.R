# # NIRMALA TIPS: 
# Very noisy - take out rare ones
# Anything above 100 ok, but to be sure go above 500 (if 100 is too noisy, switch to 500 and can't go wrong here). 
# Strength of study lies in isoform expression (we identified some rare ones), 
#       if you find reasonable consistent modules and enrichment then stay here
# Rarer you go, co-expression becomes harder because it pushes to put everything in modules

# libraries ---------------------------------------------------------------

library(tidyverse)
library(janitor)
library(RColorBrewer)
library(DESeq2)
library(WGCNA)
library(purrr)
library(patchwork)
library(ggvenn)
library(ggalluvial)
library(viridis)
library(ggpointdensity)
library(ggside)


# set theme for plots -----------------------------------------------------

theme_set(theme_bw() +
            theme(plot.title = element_text(size = 12),
                  axis.title = element_text(size = 10),
                  axis.text = element_text(size = 10),
                  strip.text = element_text(size = 10),
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 8)
                  )
          )


# data --------------------------------------------------------------------


# transcript counts
df_transcript_raw_counts <- read.csv(paste0(base_dir, "data/transcript_count_matrix.csv")) %>% 
  as_tibble() %>% 
  dplyr::rename("ensembl_transcript_id" = "X") %>% 
  rename_all(~str_remove(.x, "_MergedBam_stringtieOutput")) %>% 
  clean_names() %>% 
  rename_all(~str_remove(.x, "x")) 

# covariate data
load(paste0(base_dir, "objects/185_all_covariates_clean.RDS"))

# filter out samples that weren't in genes
df_transcript_raw_counts <- df_transcript_raw_counts %>% 
  dplyr::select(ensembl_transcript_id, contains(df_covariates_all_clean$sample))



# filter rare transcripts -------------------------------------------------

    
# ## Gandal: 0.1 TPM in at least 25% of samples
# f_tpm <- function(x) {x/sum(x) * 10^6}
# df_tpm <- df_transcript_raw_counts %>% mutate_if(is.numeric, f_tpm)

# Nirmala: 10 counts in at least 80% of samples
keep <- df_transcript_raw_counts %>% 
  mutate_if(is.numeric, ~ifelse(.x > 10, 1, 0)) %>% 
  mutate(thresh = rowSums(dplyr::select(., -ensembl_transcript_id))/(ncol(.) - 1), .before = 2) %>% 
  filter(thresh > 0.80) %>% 
  pull(ensembl_transcript_id)

# filter df_transcripts for these
df_transcript_raw_counts_filtered <- df_transcript_raw_counts %>% 
  filter(ensembl_transcript_id %in% keep) # 70k


# DESeq2 to normalize counts -------------------------
    
    
# CONVERT RAW COUNTS TO MATRIX
m_transcript_raw_counts <- df_transcript_raw_counts_filtered %>% 
  
  # convert to numeric matrix
  as.data.frame %>% 
  column_to_rownames("ensembl_transcript_id") %>% 
  as.matrix

# check that rows and columns are aligned
all(df_covariates_all_clean$sample == colnames(m_transcript_raw_counts))

# # COVARIATES: center and scale numeric variables
#     center_scale <- function(x) {(x - mean(x))/sd(x)}
#     df_covariates_scaled <- df_covariates %>% mutate_at(vars("rin_acsg", "gc_percent"), center_scale)

# CREATE DESEQ OBJECT
dds <- DESeqDataSetFromMatrix(countData = m_transcript_raw_counts,
                              colData = df_covariates_all_clean,
                              design = ~ new_dx)

# GET VSD COUNTS

# estimate DESeq size (normalization) factors
dds <- estimateSizeFactors(dds)

# extraction
m_dds_counts <- counts(dds, normalized = TRUE)

# generate vsd
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

    

# remove transcripts with low variance (relative to median expression) -----------------

# extract transformed counts from vsd object
transformed_counts <- vsd %>% assay

# pull transcript IDs
transcripts <- transformed_counts %>% rownames

# calculate median and variance for each transcript
df_transcript_variance <- map_dfr(
  .x = transcripts,
  .f = ~ tibble(
    ensembl_transcript_id = .x,
    median = median(transformed_counts[.x,]),
    sd = sd(transformed_counts[.x,]),
    n = length(transformed_counts[.x,]),
    error = qt(0.975, df = n - 1)*sd/sqrt(n),
    t_stat = t.test(transformed_counts[.x,])$statistic,
    p_val = t.test(transformed_counts[.x,])$p.value
  )
) %>% 
  mutate(ratio = sd/median) %>% 
  arrange(ratio)

# mean - variance plots
lm(sd ~ median, data = df_transcript_variance) %>% summary
df_transcript_variance %>% 
  ggplot(aes(x = median, y = sd)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Median - standard deviation relationship across transcripts")

df_transcript_variance %>% 
  ggplot(aes(x = median, y = ratio)) +
  geom_pointdensity() +
  geom_xsidehistogram() +
  geom_ysidehistogram() +
  scale_color_viridis() +
  labs(title = "Median - standard deviation:median relationship across transcripts") +
  theme_bw() 

# select threshold for removal
df_transcript_variance <- df_transcript_variance %>% 
  mutate(
    median_q1 = quantile(median, 0.25),
    ratio_q1 = quantile(ratio, 0.25),
    sd_q1 = quantile(sd, 0.25),
    
    low = ifelse(sd <= sd_q1, 1, 0)
  ) %>% 
  arrange(-low, median) 

df_transcript_variance %>% #dplyr::count(low)
  filter(low == 1) %>% 
  ggplot(aes(x = median)) +
  geom_histogram() # ok - not all low-expressed transcripts

low_var_transcripts <- df_transcript_variance %>% filter(low == 1) %>% pull(ensembl_transcript_id)

# PCA - identify covariates ---------------------------------------------------------------------

    
## RAW COUNTS

pca <- prcomp(df_transcript_raw_counts_filtered %>% 
                as.data.frame %>%
                column_to_rownames("ensembl_transcript_id"),
              center = TRUE, scale = TRUE)

df_pca <- pca$rotation %>% 
  as.data.frame %>% 
  rownames_to_column("sample") %>% 
  as_tibble() %>% 
  left_join(df_covariates_all_clean) %>% 
  mutate(library_batch = as.factor(library_batch)) %>% 
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

# plot
df_pca %>% 
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(color = "black", size = 2.75) +
  geom_point(aes(color = new_dx, text = sample), size = 2) +
  
  # label samples
  # geom_text_repel(aes(label = sample)) +
  
  # color gradient for continuous color variables
  # scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  
  # add origin lines
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  
  # give % variance in axis labels
  xlab(paste0("PC1: ", 
              df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
              "% of variance")) +
  ylab(paste0("PC2: ", 
              df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
              "% of variance")) +
  ggtitle("transcripts PCA1: no regression")

# ggplotly(p_pca, tooltip = "text")


# correlate with covariates

f_correlate_pc_covariate <- function(df_pca, df_covariates) {
  
  pcs <- df_pca %>% dplyr::select(contains("pc")) %>% colnames
  covariates <- df_covariates %>% dplyr::select(2:ncol(.)) %>% colnames
  df_pc_covar_cor <- tibble()
  
  for (pc in pcs) {
    
    for (var in covariates) {
      
      cor <- cor.test(df_pca[[pc]], df_covariates[[var]])
      estimate <- cor$estimate
      p <- cor$p.value
      
      df_tmp <- tibble(pc = pc,
                       covariate = var,
                       cor = estimate,
                       p_val = p)
      
      df_pc_covar_cor <- df_pc_covar_cor %>% bind_rows(df_tmp)
      
    }
    
  }
  
  # significant correlations
  df_pc_covar_cor %>% 
    mutate(pc = str_remove(pc, "pc") %>% as.numeric()) %>% 
    left_join(df_var %>% 
                pivot_longer(2:ncol(.), names_to = "metric", values_to = "value") %>% 
                filter(metric == "proportion_of_variance"), by = "pc") %>% 
    dplyr::select(-metric) %>% 
    relocate(value, .before = covariate) %>% 
    dplyr::rename("pc_perc_var" = "value") %>% 
    mutate(pc_perc_var = 100*pc_perc_var) %>% 
    mutate(p_adj = p.adjust(p_val, method = "BH"))
  
}

f_correlate_pc_covariate(df_pca, df_covariates_all_clean_numeric) %>% filter(p_val < 0.05 & pc_perc_var > 1)

## NORMALIZED COUNTS    

# get vsd for each gene
df_vsd <- vsd %>% assay %>% as.data.frame %>% rownames_to_column("ensembl_transcript_id") %>% as_tibble()

# run PCA on normalized counts
pca <- prcomp(df_vsd %>% 
                as.data.frame %>%
                column_to_rownames("ensembl_transcript_id"),
              center = TRUE, scale = TRUE)

df_pca <- pca$rotation %>% 
  as.data.frame %>% 
  rownames_to_column("sample") %>% 
  as_tibble() %>% 
  left_join(df_covariates_all_clean) %>% 
  mutate(library_batch = as.factor(library_batch)) %>% 
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

# plot
df_pca %>% 
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(color = "black", size = 2.75) +
  geom_point(aes(color = new_dx, text = sample), size = 2) +
  
  # label samples
  # geom_text_repel(aes(label = sample)) +
  
  # color gradient for continuous color variables
  # scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  
  # add origin lines
  #geom_vline(xintercept = 0, color = "black", lty = 2) +
  #geom_hline(yintercept = 0, color = "black", lty = 2) +
  
# give % variance in axis labels
xlab(paste0("PC1: ", 
            df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
            "% of variance")) +
  ylab(paste0("PC2: ", 
              df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
              "% of variance")) +
  ggtitle("transcripts PCA1: normalized counts + no regression")

# correlate
f_correlate_pc_covariate <- function(df_pca, df_covariates) {
  
  pcs <- df_pca %>% dplyr::select(contains("pc")) %>% colnames
  covariates <- df_covariates %>% dplyr::select(2:ncol(.)) %>% colnames
  df_pc_covar_cor <- tibble()
  
  for (pc in pcs) {
    
    for (var in covariates) {
      
      cor <- cor.test(df_pca[[pc]], df_covariates[[var]])
      estimate <- cor$estimate
      p <- cor$p.value
      
      df_tmp <- tibble(pc = pc,
                       covariate = var,
                       cor = estimate,
                       p_val = p)
      
      df_pc_covar_cor <- df_pc_covar_cor %>% bind_rows(df_tmp)
      
    }
    
  }
  
  # significant correlations
  df_pc_covar_cor %>% 
    mutate(pc = str_remove(pc, "pc") %>% as.numeric()) %>% 
    left_join(df_var %>% 
                pivot_longer(2:ncol(.), names_to = "metric", values_to = "value") %>% 
                filter(metric == "proportion_of_variance"), by = "pc") %>% 
    dplyr::select(-metric) %>% 
    relocate(value, .before = covariate) %>% 
    dplyr::rename("pc_perc_var" = "value") %>% 
    mutate(pc_perc_var = 100*pc_perc_var) %>% 
    mutate(p_adj = p.adjust(p_val, method = "BH"))
  
}
f_correlate_pc_covariate(df_pca, df_covariates_all_clean_numeric) %>% filter(p_adj < 0.05 & pc_perc_var > 5)



# regress significant covariates ------------------------------------------------------

    
# REGRESS COVARIATES

df_vsd_t <- vsd %>% assay %>% t %>% as.data.frame %>% as_tibble

# fit lm to gene using covariates (GCpercent, RIN, Race), get residuals for all genes
df_vsd_resids_transcripts <- tibble()

for (i in 1:ncol(df_vsd_t)) {
  
  fit <- lm(df_vsd_t[[i]] ~ 
              df_covariates_all_clean$rin_acsg +
              df_covariates_all_clean$mapped_percent +
              df_covariates_all_clean$g_cpercent +
              df_covariates_all_clean$race) %>% summary
  
  df_tmp <- residuals(fit) %>% # extract residuals
    as.data.frame %>% 
    t %>% 
    as_tibble %>% 
    mutate(ensembl_transcript_id = df_vsd_t[i] %>% colnames, .before = 1) # add gene ID for rows
  
  df_vsd_resids_transcripts <- df_vsd_resids_transcripts %>% bind_rows(df_tmp)
  
}

df_vsd_resids_transcripts <- df_vsd_resids_transcripts  %>% rename_at(2:ncol(.), ~df_covariates_all_clean$sample) # add column names as samples

# write to csv
write.csv(df_vsd_resids_transcripts, 
          file = "outputs/20221122_185samples_70ktranscripts_vsd_rin_mp_gc_race_resids.csv",
          row.names = FALSE)


# TRANSPOSE expression data for further analysis, rename columns & add ID row

df_vsd_resids_transcripts_t <- df_vsd_resids_transcripts %>%
  
  # transpose & convert to tibble
  dplyr::select(-1) %>%
  t %>%
  as.data.frame %>%
  as_tibble %>%
  
  # rename columns with gene id from original df
  rename_all(~ df_vsd_resids_transcripts %>% pull(ensembl_transcript_id)) %>%
  
  # add sample id col using colnames of original df (to maintain order)
  dplyr::mutate(sample = df_vsd_resids_transcripts %>% dplyr::select(-1) %>% colnames, .before = 1)

# save object to move to cluster
save(df_vsd_resids_transcripts_t, 
     file = "objects/20221122_185samples_70ktranscripts_vsd_rin_mp_gc_race_resids_t.Rdata")


# PCA on regressed data ---------------------------------------------------


      # run PCA
        pca <- prcomp(df_vsd_resids_transcripts %>% 
                        as.data.frame %>%
                        column_to_rownames("ensembl_transcript_id"),
                      center = TRUE, scale = TRUE)
        
        df_pca <- pca$rotation %>% 
          as.data.frame %>% 
          rownames_to_column("sample") %>% 
          as_tibble() %>% 
          left_join(df_covariates_all_clean) %>% 
          mutate(library_batch = as.factor(library_batch)) %>% 
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
        
      # plot pca
        p_pca <- df_pca %>% 
          ggplot(aes(x = pc1, y = pc2)) +
          geom_point(color = "black", size = 2.75) +
          geom_point(aes(color = new_dx, text = sample), size = 2) +
          
          # label samples
          # geom_text_repel(aes(label = sample)) +
          
          # color gradient for continuous color variables
          # scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
          
          # add origin lines
          geom_vline(xintercept = 0, color = "black", lty = 2) +
          geom_hline(yintercept = 0, color = "black", lty = 2) +
          
          # give % variance in axis labels
          xlab(paste0("PC1: ", 
                      df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                      "% of variance")) +
          ylab(paste0("PC2: ", 
                      df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                      "% of variance")) +
          ggtitle("PCA on transformed, residualized transcript counts")
        
      # plot eigenvalues
        p_var <- df_var %>% 
          ggplot(aes(x = pc, y = proportion_of_variance)) +
          geom_bar(stat = "identity") +
          geom_hline(aes(yintercept = 0.01), lty = 2, color = "black")
        
        p_pca + p_var
        
      # correlate
        df_pc_covar_cor <- f_correlate_pc_covariate(df_pca, df_covariates_all_clean_numeric)
        
      # plot correlations
        df_pc_covar_cor %>% 
          filter(pc_perc_var > 1) %>% 
          mutate(p_adj = ifelse(p_adj == 1, NA, p_adj)) %>% 
          mutate(labels = ifelse(p_adj < 0.05, formatC(p_adj, format = "e", digits = 2), "")) %>% 
          
          ggplot(aes(x = factor(pc), y = covariate)) +
          geom_tile(aes(fill = -log10(p_adj))) +
          geom_text(aes(label = labels), color = "white", size = 3) +
          geom_vline(aes(xintercept = 1.5), lty = 2, color = "black") +
          scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
          labs(x = "PC (> 1% explained var)") +
          theme(axis.text = element_text(size = 15),
                axis.title = element_text(size = 20),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 15),
                legend.position = "bottom") +
          ggtitle("correlate PCs with covariates")
        

# transcript PCA ----------------------------------------------------------


        pca_transcripts <- prcomp(df_vsd_resids_transcripts_t %>% as.data.frame %>% column_to_rownames("sample"))
        
        df_pca_transcripts <- pca_transcripts$rotation %>% 
          as.data.frame %>% 
          rownames_to_column("ensembl_transcript_id") %>% 
          as_tibble() %>% 
          clean_names()
        
        df_var_transcripts <- summary(pca_transcripts) %>% 
          .$importance %>% 
          as.data.frame %>% 
          rownames_to_column("metric") %>% 
          as_tibble() %>% 
          pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
          mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
          pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
          clean_names()
        
        df_pca_transcripts %>% 
          ggplot(aes(x = pc1, y = pc2)) +
          geom_point(color = "black", size = 2.75, shape = 1) +
          
          # label samples
          # geom_text_repel(aes(label = sample)) +
          
          # add origin lines
          geom_vline(xintercept = 0, color = "black", lty = 2) +
          geom_hline(yintercept = 0, color = "black", lty = 2) +
          
          # give % variance in axis labels
          xlab(paste0("PC1: ", 
                      df_var_transcripts %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                      "% of variance")) +
          ylab(paste0("PC2: ", 
                      df_var_transcripts %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                      "% of variance"))
        
    

# select soft-thresholding power ------------------------------------------

    
    # RUN SFTS
        
        # enable parallel processing
        disableWGCNAThreads()
        doParallel::registerDoParallel()
        
        # use matrix for these steps
        m_vsd_resids_transcript_t <- df_vsd_resids_transcripts_t %>% as.data.frame %>% column_to_rownames("sample")
        
        # choose a set of soft-thresholding powers
        powers <- 1:16
        
        # Call the network topology analysis function
        sft <- pickSoftThreshold(m_vsd_resids_transcript_t, powerVector = powers, verbose = 5)
        
        # plot the results
        df_sft <- sft$fitIndices %>% as_tibble %>% clean_names
        
        # save object
        save(df_sft, file = "objects/20221122_185samples_70ktranscripts_vsd_mp_rin_gc_rac_resids_sft.RDS")
        
        # scale-free topology fit index as a function of the soft-thresholding power
        p1 <- df_sft %>% 
          ggplot(aes(x = power, y = -sign(slope)*sft_r_sq)) +
          geom_text(aes(label = power), size = 5) +
          
          # geom_hline(yintercept = 0.8, lty = 2, color = "red") +
          geom_hline(yintercept = 0.9, lty = 2, color = "red") +
          ylim(c(0, 1)) +
          xlab("Soft threshold (power)") +
          ylab("Scale free topology model fit (signed R^2)") 
    
    # PLOT
        
        # Median connectivity as a function of the soft-thresholding power
        p2 <- df_sft %>% 
          ggplot(aes(x = power, y = median_k)) +
          geom_text(aes(label = power), size = 5) +
          geom_hline(yintercept = 100, lty = 2, color = "black") +
          
          xlab("Soft threshold (power)") +
          ylab("Median connectivity") # check that mean connectivity remains reasonably high (in the 100s) or above
        
        p1 + p2
        
        # use WGCNA suggestion
        soft_power <- sft$powerEstimate

        
        
# adjacency matrix, TOM, clustering ------------------------------------
        
        doParallel::registerDoParallel()
        
    # load & configure data (if starting here)
        
        # load
        #load("objects/20221108_185samples_18ktranscripts_vsd_rin_mp_gc_race_resids_t.Rdata")
        load("objects/20221102_185samples_44ktranscripts_vsd_rin_mp_gc_race_resids_t.Rdata")
        
        # transpose to matrix
        m_vsd_resids_transcripts_t <- df_vsd_resids_transcripts_t %>% as.data.frame %>% column_to_rownames("sample")
        
        rownames(m_vsd_resids_transcripts_t) <- df_vsd_resids_transcripts_t %>% pull(sample)
        colnames(m_vsd_resids_transcripts_t) <- df_vsd_resids_transcripts_t %>% dplyr::select(-sample) %>% colnames()
        
        # select sft here
        soft_power <- 3
        
        # co-expression similarity and adjacency
        adjacency <- adjacency(m_vsd_resids_transcripts_t, power = soft_power)
        
        # topological overlap matrix (TOM)
        TOM <- TOMsimilarity(adjacency)
        diss_TOM <- 1 - TOM
        
        # clustering using TOM
        
        # call the hierarchical clustering function
        gene_tree <- hclust(as.dist(diss_TOM))
        
        # plot the resulting clustering tree
        # par(mfrow = c(1, 1))
        # plot(gene_tree, labels = FALSE, xlab = "")
        
    
        
# module assignment -------------------------------------------------------
        
        
        doParallel::registerDoParallel()
        
    # load data
        load("objects/20221108_185samples_44ktranscripts_vsd_rin_mp_gc_race_resids_t.Rdata")
        
    # LOOP to assign modules based on different MIN_MOD_SIZE and CUT_HEIGHT
        
        min_mod_size <- c(30, 50)
        cut_heights <- c(0.15, 0.20, 0.25)
        df_modules_transcripts <- tibble()
        
        for (i in min_mod_size) {
          
          print(paste0("assigning modules for minimum module size = ", i))
          
          dynamic_mods <- cutreeDynamic(dendro = gene_tree,
                                        distM = diss_TOM,
                                        # deepSplit = split, # use default deepSplit param
                                        pamRespectsDendro = FALSE,
                                        minClusterSize = i)
          dynamic_colors <- labels2colors(dynamic_mods)
          
          for (j in cut_heights) {
            
            print(paste0("merging modules at cut_height = ", j))
            
            # call an automatic merging function using cutHeight
            merge <- mergeCloseModules(m_vsd_resids_transcripts_t,
                                       dynamic_mods,
                                       cutHeight = j,
                                       verbose = 3)
            merged_mods <- merge$colors
            merged_colors <- labels2colors(merged_mods)
            
            # create df
            df_tmp <- tibble(min_size = i,
                             cut_height = j,
                             ensembl_transcript_id = df_vsd_resids_transcripts$ensembl_transcript_id,
                             dynamic_mod_color = dynamic_colors,
                             merged_mod_color = merged_colors)
            
            df_modules_transcripts <- df_modules_transcripts %>% bind_rows(df_tmp)
            
          }
          
          
        }
        
        df_modules_transcripts %>% dplyr::count(min_size, cut_height)
        save(df_modules_transcripts, 
             file = "objects/20221108_185samples_44ktranscripts_vsd_rin_mp_gc_race_resids_sft3_minSize_cutHeight_mods.RDS")
        
        
        
        
# combine and plot module sizes -------------------------------------------
        
    
# LOAD DATA FROM BIOWULF            
        load("objects/20221123_185samples_70ktranscripts_vsd_rin_mp_gc_race_SFT8_minSize_cutHeight_mods.RDS")
        df_modules_transcripts <- df_modules_transcripts %>% 
          dplyr::mutate(sft = 8, .before = 1) %>% 
          mutate(mod_set = paste0("sft", sft, "_minSize", min_size, "_cutHeight", cut_height), .before = 1)

# RENUMBER
        df_dynamic_mod_no <- df_modules_transcripts %>% 
          group_by(sft, min_size, cut_height) %>% 
          dplyr::count(dynamic_mod_color) %>% 
          arrange(sft, min_size, cut_height, desc(dynamic_mod_color == "grey"), desc(n)) %>% 
          mutate(dynamic_module = 0:(n() - 1) %>% factor(levels = 0:(n() - 1)))
        
        df_merged_mod_no <- df_modules_transcripts %>% 
          group_by(sft, min_size, cut_height) %>% 
          dplyr::count(merged_mod_color) %>% 
          arrange(sft, min_size, cut_height, desc(merged_mod_color == "grey"), desc(n)) %>% 
          mutate(merged_module = 0:(n() - 1) %>% factor(levels = 0:(n() - 1)))
        
        df_modules_transcripts <- df_modules_transcripts %>% 
          left_join(df_dynamic_mod_no %>% dplyr::select(-n)) %>% 
          left_join(df_merged_mod_no %>% dplyr::select(-n)) %>% 
          dplyr::rename("module" = "merged_module")
        
# COMBINE WITH GENES
        load("objects/transcript_gene_go_term.RDS")
        
        df_modules_transcripts <- df_modules_transcripts %>% 
          left_join(df_transcript_gene_go_term %>% 
                      dplyr::select(ensembl_gene_id, ensembl_transcript_id) %>% 
                      distinct(),
                    by = "ensembl_transcript_id")
        
        # save object
        save(df_modules_transcripts, 
             file = "objects/20221123_185samples_70ktranscripts_vsd_rin_mp_gc_race_SFT8_minSize_cutHeight_mods.RDS")
        
        
  ## PLOT
        
        # load("objects/20221108_185samples_18ktranscripts_vsd_rin_mp_gc_rac_sft4_minSize_cutHeight_mods.RDS")
        
        # transcript level
        df_module_sizes <- df_modules_transcripts %>%
          group_by(sft, min_size, cut_height) %>%
          dplyr::count(module)
        
        df_module_sizes %>%
          # mutate(label = ifelse(as.numeric(as.character(module)) %% 5 == 0, module, "")) %>% 
          ggplot(aes(x = module, y = n)) +
          geom_point(shape = 1, size = 2) +
          facet_grid(min_size ~ cut_height, scales = "free") +
          ggtitle("sgACC 70k transcripts, sft8 module assignments") +
          theme(axis.text.x = element_text(angle = 90))
        
        # gene level
        df_modules_transcripts %>% pull(ensembl_gene_id) %>% unique %>% length()
        
        df_module_sizes <- df_modules_transcripts %>%
          group_by(sft, min_size, cut_height) %>%
          dplyr::select(ensembl_gene_id, module) %>% 
          distinct() %>% 
          dplyr::count(module)
        
        df_module_sizes %>%
          ggplot(aes(x = module, y = n)) +
          geom_point(shape = 1, size = 2) +
          facet_grid(min_size ~ cut_height, scales = "free") +
          ggtitle("sgACC 70k transcripts --> 21k genes, sft8 module assignments") +
          theme(axis.text.x = element_text(angle = 90))
        


# compare transcript level to gene level --> any new info?? ---------------


## LOAD DATA
        
        load("objects/20221123_185samples_70ktranscripts_vsd_rin_mp_gc_race_SFT8_minSize_cutHeight_mods.RDS") # df_modules_transcripts
        df_modules_transcripts <- df_modules_transcripts %>% filter(mod_set == "sft8_minSize30_cutHeight0.2")
        
        load("objects/20220613_185samples_14kgenes_vsd_RIN_mappedPerc_GC_race_resids_sft10_11_minSize_cutHeight_mods.RDS") # df_modules
        df_modules <- df_modules %>% filter(mod_set == "sft10_minSize30_cutHeight0.15")

        
        
## FIND GREY MOD OVERLAP
        
        transcript_grey_mod_genes <- df_modules_transcripts %>% 
          filter(module == 0 & !is.na(ensembl_gene_id)) %>% 
          pull(ensembl_gene_id) %>% unique
        
        gene_grey_mod_genes <- df_modules %>% 
          filter(merged_module == 0 & !is.na(ensembl_gene_id)) %>% 
          pull(ensembl_gene_id) %>% unique
        
        length(transcript_grey_mod_genes)
        length(gene_grey_mod_genes)
        
        length(intersect(transcript_grey_mod_genes, gene_grey_mod_genes))
        
        # venn diagram
        p_grey <- ggvenn(list(transcript_level = transcript_grey_mod_genes, gene_level = gene_grey_mod_genes)) +
          ggtitle("Unassigned (grey module) genes at transcript and gene level")
        
## FIND ASSIGNED GENE OVERLAP
        
        transcript_assigned_genes <- df_modules_transcripts %>% 
          filter(module != 0 & !is.na(ensembl_gene_id)) %>% 
          pull(ensembl_gene_id) %>% unique
        
        gene_assigned_genes <- df_modules %>% 
          filter(merged_module != 0 & !is.na(ensembl_gene_id)) %>% 
          pull(ensembl_gene_id) %>% unique
        
        length(transcript_assigned_genes)
        length(gene_assigned_genes)
        
        length(intersect(transcript_assigned_genes, gene_assigned_genes))
        
        # venn diagram
        p_assigned <- ggvenn(list(transcript_level = transcript_assigned_genes, gene_level = gene_assigned_genes)) +
          ggtitle("Assigned genes at transcript and gene level")
        
        p_grey / p_assigned
        
# WHERE DID RARE TRANSCRIPTS END UP??
        
        # define rare - start with mean counts < 50 across all samples
        rare_transcripts <- df_transcript_raw_counts_filtered %>% 
          pivot_longer(2:ncol(.), names_to = "sample", values_to = "counts") %>% 
          group_by(ensembl_transcript_id) %>% 
          summarise(mean = mean(counts)) %>% 
          filter(mean < 100) %>% 
          pull(ensembl_transcript_id)
          
        df_modules_transcripts %>% 
          filter(ensembl_transcript_id %in% rare_transcripts) %>% 
          dplyr::count(module)
      
      
# ALLUVIAL
        
        load("objects/transcript_gene_go_term.RDS")
        
        # get alluvial format
        df_alluvial <- df_modules %>% 
          left_join(df_transcript_gene_go_term %>% 
                      dplyr::select(ensembl_gene_id, ensembl_transcript_id) %>% 
                      distinct(),
                    by = "ensembl_gene_id") %>%  # 1164 genes without transcripts
          dplyr::select(mod_set, sft, min_size, cut_height, ensembl_gene_id, ensembl_transcript_id, merged_module) %>% 
          dplyr::rename("module" = "merged_module") %>% 
          mutate(mod_set = paste0("genes_", mod_set)) %>% 
            
          bind_rows(df_modules_transcripts %>% 
                      select(mod_set, sft, min_size, cut_height, ensembl_gene_id, ensembl_transcript_id, module) %>% 
                      mutate(mod_set = paste0("transcripts_", mod_set))) %>% 
          mutate(mod_set = factor(mod_set, levels = unique(.$mod_set))) %>% 
          distinct()

        # rename and add colors
        qual_col_pals <- brewer.pal.info %>% 
          rownames_to_column("palette") %>% 
          filter(category == "qual", !grepl("Pastel", palette))
        
        col_vector <- c("#636363", unlist(mapply(brewer.pal, 
                                                 (qual_col_pals$maxcolors - 1), 
                                                 qual_col_pals$palette)))
        df_colors <- tibble(module = levels(df_alluvial$module),
                            color = col_vector[1:length(levels(df_alluvial$module))]) %>% 
          mutate(module = factor(module, levels = 0:max(as.numeric(module))))
        
        df_alluvial <- df_alluvial %>% 
          dplyr::select(ensembl_transcript_id, mod_set, module) %>%
          left_join(df_colors) %>% 
          distinct()

        # plot
        df_alluvial %>% 
          
          # mutate(mod_set = factor(mod_set, levels = c("transcripts_sft8_minSize30_cutHeight0.2", "genes_sft10_minSize30_cutHeight0.15"))) %>% 
          
          ggplot(aes(x = mod_set, stratum = module, alluvium = ensembl_transcript_id)) +
          stat_stratum() +
          geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
          stat_flow(aes(fill = I(color))) +
          ggtitle("70k transcripts")
  
# TRANSCRIPTS FROM SAME GENE IN MORE THAN ONE MODULE?
        
        multi_mod_genes <- df_modules_transcripts %>% 
          dplyr::count(ensembl_gene_id, module) %>% 
          filter(n > 1 & !is.na(ensembl_gene_id)) %>% 
          arrange(-n) %>% 
          pull(ensembl_gene_id)
              
        df_modules_transcripts %>% 
          filter(mod_set == "sft10_minSize30_cutHeight0.1" & ensembl_gene_id == "ENSG00000169398") %>%
          dplyr::count(module)
          
        df_alluvial %>% 
          filter(ensembl_gene_id %in% multi_mod_genes[1:100]) %>% 
          
          #filter(mod_set %in% c("transcripts_sft10_minSize30_cutHeight0.1", "genes_sft10_minSize30_cutHeight0.15")) %>% 
          left_join(df_colors) %>% 
          mutate(mod_set = factor(ifelse(grepl("transcript", mod_set), "transcripts", "gene"),  
                                  levels = c("transcripts", "gene")
                                  )
                 ) %>% 
          
          ggplot(aes(x = mod_set, stratum = module, alluvium = ensembl_transcript_id)) +
          stat_stratum() +
          geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
          stat_flow(aes(fill = I(color))) +
          facet_wrap(~ensembl_gene_id, scales = "free") +
          ggtitle("TPM > 1; 44k transcripts") +
          theme(strip.text = element_text(size = 8))
        
        
        
# relate modules to external clinical traits (basic correlation) ------------------------------
        
        
    # QUANTIFY MODULE-TRAIT ASSOCIATIONS
        
        # load data 
        load("objects/20221123_185samples_70ktranscripts_vsd_rin_mp_gc_race_SFT8_minSize_cutHeight_mods.RDS") # df_modules_transcripts
        df_modules_transcripts <- df_modules_transcripts %>% filter(mod_set == "sft8_minSize30_cutHeight0.2")
        load("objects/20221122_185samples_70ktranscripts_vsd_rin_mp_gc_race_resids_t.Rdata")
        
        # define numbers of genes and samples
        m_vsd_resids_transcripts_t <- df_vsd_resids_transcripts_t %>%
          as.data.frame %>% 
          column_to_rownames("sample") %>% 
          as.matrix
        
        n_transcripts <- ncol(m_vsd_resids_transcripts_t)
        n_samples <- nrow(m_vsd_resids_transcripts_t)
        
        # get numeric covariates (df_covariates_numeric loaded in `data`)
        numeric_covariates <- df_covariates_all_clean_numeric %>% column_to_rownames("sample")
        
        # recalculate MEs with their color labels
        merged_colors <- df_modules_transcripts$merged_mod_color
        MEs0 <- moduleEigengenes(m_vsd_resids_transcripts_t, merged_colors)$eigengenes
        MEs <- orderMEs(MEs0)
        module_trait_cor <- cor(MEs, numeric_covariates, use = "p")
        module_trait_p_value <- corPvalueStudent(module_trait_cor, n_samples)
        
        # create df of module-covariate correlations and associated p-values
        sum(rownames(module_trait_cor) != rownames(module_trait_p_value)) # check that modules are in same order
        
        # plot
        source("functions/f_plot_module_trait_correlation.R")
        
        f_mod_trait_cor(m_sample_exp = m_vsd_resids_transcripts_t,
                        numeric_covariates = numeric_covariates,
                        merged_colors = merged_colors) +
          theme(axis.text.x = element_text(angle = 45, hjust = 0.9)) +
          ggtitle("No significant module-trait correlations")
        
        

# isoform functional enrichment -------------------------------------------

        
        # load transcript-GO data
        load("objects/transcript_gene_go_term.RDS")
        
        # source function
        source("functions/f_top_GO_modules.R")
        load("objects/df_go_full_definitions.Rdata") # full go definitions
        
        # load in module data

        # map transcript ID to gene ID
        df_modules_genelevel <- df_modules_transcripts %>% 
          dplyr::select(-ensembl_transcript_id) %>% 
          distinct()
      
    # PULL ALL WGCNA GENES AFTER FILTERING FOR BACKGROUND
        l_gene_universe <- df_modules_genelevel$ensembl_gene_id %>% unique
        
    # CREATE GO ANNOTATION
        
        # get gene info from ensembl database
        hsapiens_ensembl <- useEnsembl(biomart = "genes",
                                       dataset = "hsapiens_gene_ensembl")
        
        # map each gene to GO term
        m_go_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "namespace_1003"),
                          filters = "ensembl_gene_id",
                          values = l_gene_universe,
                          mart = hsapiens_ensembl)
        
        # build annotation list
        l_gene_2_GO <- m_go_ids[,c(1,3)] %>% unstack
        
        
    # RUN ENRICHMENT FUNCTION ON ALL MODULES IN ALL MODULE SETS
        doParallel::registerDoParallel()
        
        mod_sets <- df_modules_genelevel$mod_set %>% unique
        df_mods_go <- tibble()
        for (set in mod_sets) {
          
          print(paste0("determining functional enrichment for module set ", set))
          
          df_mod_set <- df_modules_genelevel %>% filter(mod_set == set)
          modules <- df_mod_set %>% filter(!is.na(module)) %>% arrange(module) %>% pull(module) %>% unique
          n_mods <- length(modules) - 1
          
          df_tmp <- 0:n_mods %>%
            map_dfr(~f_top_GO_modules(l_gene_module = df_mod_set %>%
                                        filter(module == .x) %>%
                                        pull(ensembl_gene_id),
                                      l_gene_universe = l_gene_universe,
                                      m_go_ids = m_go_ids,
                                      l_gene_2_GO = l_gene_2_GO,
                                      go_ontology = "BP") %>%
                      mutate(module = .x, .before = 1)
            ) %>%
            clean_names() %>%
            dplyr::select(-term) %>%
            left_join(df_go_defs, by = "go_id") %>%
            dplyr::select(module, go_id, term, everything()) %>%
            distinct() %>%
            mutate(mod_set = set, .before = 1)
          
          df_mods_go <- df_mods_go %>% bind_rows(df_tmp)
          
        }
        
        
        # EXPORT TO EXCEL & SAVE OBJECT
        write.csv(df_mods_go,
                  file = "outputs/20221127_185samples_70ktranscripts_vsd_rin_mp_gc_rac_sft8_GO_RES.csv")
        
        save(df_mods_go, file = "objects/20221127_185samples_70ktranscripts_vsd_rin_mp_gc_rac_sft8_GO_RES.RDS")

        
    # PLOT
        require(tidytext)
        df_mods_go %>% 
          
          group_by(module) %>% 
          arrange(module, p_adj) %>% 
          dplyr::top_n(n = 5) %>% 
          
          ggplot(aes(x = -log10(p_adj), y = reorder_within(str_wrap(term, width = 30), -p_adj, module))) +
          geom_point(aes(size = significant/annotated, color = -log10(p_adj))) +
          labs(y = "") +
          facet_wrap(~module, scales = "free") +
          scale_y_discrete(labels = function(x) gsub("\\_.*", "", x)) +
          
          geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
          geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
          geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +
          scale_color_gradientn(colors = brewer.pal(9, "Reds")[3:9]) +
          theme(plot.title = element_text(size = 18),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 8),
                strip.text = element_text(size = 10),
                legend.title = element_text(size = 15),
                legend.text = element_text(size = 12)) +
          ggtitle("70k transcripts; sft8_minSize30_cutHeight0.2")
        
        

# module eigengene-diagnosis correlation ----------------------------------


# PULL MERGED MODULE EIGENGENES
        
        merged_colors <- df_modules_transcripts %>% 
          dplyr::select(module, ensembl_transcript_id) %>% 
          distinct() %>% 
          pull(module)
        
        # assign colors to numbers
        df_color_mod <- df_modules_transcripts %>% 
          dplyr::select(module, merged_mod_color, ensembl_transcript_id) %>% 
          distinct() %>% 
          dplyr::count(merged_mod_color) %>%
          arrange(desc(merged_mod_color == "grey"), desc(n)) %>%
          dplyr::rename("color" = "merged_mod_color") %>%
          mutate(module = 0:(n() - 1) %>% factor(levels = 0:(n() - 1)), .before = 1)
        
        
        
  # GET EIGENGENES
        
        MEs <- moduleEigengenes(m_vsd_resids_transcripts_t, merged_colors[, drop = TRUE])$eigengenes
        
        df_MEs <- MEs %>%
          as.data.frame %>%
          rename_all(~str_remove(.x, "ME")) %>%
          rownames_to_column("sample") %>%
          as_tibble()
        
        df_MEs_long <- df_MEs %>%
          pivot_longer(2:ncol(.), names_to = "module", values_to = "kme") %>%
          left_join(df_color_mod, by = "module") %>%
          left_join(df_covariates_all_clean, by = "sample") %>% 
          mutate(new_dx = factor(new_dx, levels = c("control", "bipolar", "mdd", "schizo"))) %>% 
          mutate(module = factor(module, levels = 0:max(as.numeric(module))))
        
  # PLOT EIGENGENES ACROSS DX
        df_MEs_long %>%
          ggplot(aes(x = module, y = kme)) +
          geom_boxplot(aes(color = new_dx))
        
        df_MEs_long %>%
          ggplot(aes(x = kme, fill = new_dx)) +
          geom_density(alpha = 0.5) +
          facet_wrap(~module) +
          theme(plot.title = element_text(size = 18),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 12),
                strip.text = element_text(size = 12),
                legend.title = element_text(size = 15),
                legend.text = element_text(size = 12)) +
          ggtitle("module eigengene distributions by diagnosis in sft8_minSize30_cutHeight0.2")
        
# KS TEST FOR DIFFERENCES IN DISTRIBUTIONS
        
        df_ks_dx <- tibble()
        
        for (i in 1:levels(df_MEs_long$module)) {
          
          kme_control <- df_MEs_long %>% filter(module == i & new_dx == "control") %>% pull(kme)
          kme_bipolar <- df_MEs_long %>% filter(module == i & new_dx == "bipolar") %>% pull(kme)
          kme_mdd <- df_MEs_long %>% filter(module == i & new_dx == "mdd") %>% pull(kme)
          kme_schizo <- df_MEs_long %>% filter(module == i & new_dx == "schizo") %>% pull(kme)
          
          ks_bipolar <- ks.test(kme_control, kme_bipolar)
          ks_mdd <- ks.test(kme_control, kme_mdd)
          ks_schizo <- ks.test(kme_control, kme_schizo)
          
          df_tmp <- tibble(module = i,
                           ks_bipolar = ks_bipolar$p.value,
                           ks_mdd = ks_mdd$p.value,
                           ks_schizo = ks_schizo$p.value)
          
          df_ks_dx <- df_ks_dx %>% bind_rows(df_tmp)
          
        }
        
        
        df_ks_dx %>% print(n = nrow(.))
 
# KRUSKAL WALLIS --> PAIRWISE WILCOX
        
        df_MEs_long %>% 
          group_by(module) %>% 
          nest() %>% 
          
          # kruskal-wallis
          mutate(kw = map(.x = data,
                          .f = ~ kruskal.test(kme ~ new_dx, data = .x))) %>% 
          mutate(kw_p_val = map(.x = kw,
                                .f = ~ .x$p.value)) %>% 
          unnest(cols = c(kw_p_val))
        
# LINEAR REGRESSION OF MES AND DIAGNOSIS
        
        df_lm_mods <- tibble()
        for (i in levels(df_MEs_long$module)) {
          
          fit <- lm(kme ~ new_dx, data = df_MEs_long %>% filter(module == i))
          sum <- fit %>% summary()
          
          df_tmp <- fit %>%
            tidy %>%
            filter(term != "(Intercept)") %>%
            mutate(term = str_remove(term, "new_dx")) %>%
            dplyr::rename("term_p_val" = "p.value", "t_stat" = "statistic", "coefficient" = "estimate") %>%
            clean_names %>%
            mutate(module = i, .before = 1) %>%
            mutate(adj_r_sq = sum$adj.r.squared,
                   model_p_val =  pf(sum$fstatistic[1], sum$fstatistic[2], sum$fstatistic[3], lower.tail = FALSE),
                   .before = 2)
          
          df_lm_mods <- df_lm_mods %>% bind_rows(df_tmp)
          
        }
        
        # combine with model on all modules combined
        fit_all <- lm(kme ~ new_dx, data = df_MEs_long)
        sum_all <- fit_all %>% summary()
        
        df_lm_mods <- fit_all %>%
          tidy %>%
          filter(term != "(Intercept)") %>%
          mutate(term = str_remove(term, "new_dx")) %>%
          dplyr::rename("term_p_val" = "p.value", "t_stat" = "statistic", "coefficient" = "estimate") %>%
          clean_names %>%
          mutate(module = "all", .before = 1) %>%
          mutate(adj_r_sq = sum$adj.r.squared,
                 model_p_val =  pf(sum_all$fstatistic[1], sum_all$fstatistic[2], sum_all$fstatistic[3], lower.tail = FALSE),
                 .before = 2) %>%
          bind_rows(df_lm_mods) %>% 
          mutate(p_adj = p.adjust(term_p_val, method = "BH")) %>% 
          distinct
        
        df_lm_mods %>% 
          filter(term_p_val < 0.05) %>% 
          print(n = nrow(.))
        
        
# LINEAR REGRESSION OF MES AND DIAGNOSIS WITH COVARIATES?
        
        df_lm_mods <- tibble()
        for (i in levels(df_MEs_long$module)) {
          
          fit <- lm(kme ~ new_dx*gender, data = df_MEs_long %>% filter(module == i))
          sum <- fit %>% summary()
          
          df_tmp <- fit %>%
            tidy %>%
            filter(term != "(Intercept)") %>%
            mutate(term = str_remove(term, "new_dx")) %>%
            dplyr::rename("term_p_val" = "p.value", "t_stat" = "statistic", "coefficient" = "estimate") %>%
            clean_names %>%
            mutate(module = i, .before = 1) %>%
            mutate(adj_r_sq = sum$adj.r.squared,
                   model_p_val =  pf(sum$fstatistic[1], sum$fstatistic[2], sum$fstatistic[3], lower.tail = FALSE),
                   .before = 2)
          
          df_lm_mods <- df_lm_mods %>% bind_rows(df_tmp)
          
        }
        
        # combine with model on all modules combined
        fit_all <- lm(kme ~ new_dx*gender, data = df_MEs_long)
        sum_all <- fit_all %>% summary()
        
        df_lm_mods <- fit_all %>%
          tidy %>%
          filter(term != "(Intercept)") %>%
          mutate(term = str_remove(term, "new_dx")) %>%
          dplyr::rename("term_p_val" = "p.value", "t_stat" = "statistic", "coefficient" = "estimate") %>%
          clean_names %>%
          mutate(module = "all", .before = 1) %>%
          mutate(adj_r_sq = sum$adj.r.squared,
                 model_p_val =  pf(sum_all$fstatistic[1], sum_all$fstatistic[2], sum_all$fstatistic[3], lower.tail = FALSE),
                 .before = 2) %>%
          bind_rows(df_lm_mods) %>% 
          mutate(p_adj = p.adjust(term_p_val, method = "BH")) %>% 
          distinct
        
        df_lm_mods %>% 
          filter(term_p_val < 0.05) %>% 
          print(n = nrow(.))
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    