
######################################################################################

# Identify qSVs using transcript counts and qsvaR to regress in step 01

######################################################################################


# libraries ---------------------------------------------------------------


  library(tidyverse)
  library(janitor)
  library(qsvaR)
  library(biomaRt)
  library(RColorBrewer)
  library(patchwork)
  library(WGCNA)
  library(DESeq2)



# transcript data ---------------------------------------------------------


# TRANSCRIPT COUNTS
  df_transcript_raw_counts <- read.csv("data/transcripts/transcript_count_matrix.csv") %>% 
    as_tibble() %>% 
    dplyr::rename("ensembl_transcript_id" = "X") %>% 
    rename_all(~str_remove(.x, "_MergedBam_stringtieOutput")) %>% 
    clean_names() %>% 
    rename_all(~str_remove(.x, "x")) 

# COVARIATES
  load("data/covariates/185_all_covariates_clean.Rdata")
  m_metadata <- df_covariates_all_clean_numeric %>% 
    as.data.frame %>% 
    column_to_rownames("sample") %>% 
    as.matrix

# FILTER COUNTS FOR ONLY SAMPLES THAT ARE IN GENE DATA
  df_transcript_raw_counts <- df_transcript_raw_counts %>% 
    dplyr::select(ensembl_transcript_id, contains(df_covariates_all_clean$sample))
  
  

# remove low count transcripts --------------------------------------------


# Nirmala: 10 counts in at least 80% of samples
  keep <- df_transcript_raw_counts %>% 
    mutate_if(is.numeric, ~ifelse(.x > 10, 1, 0)) %>% 
    mutate(thresh = rowSums(dplyr::select(., -ensembl_transcript_id))/(ncol(.) - 1), .before = 2) %>% 
    filter(thresh > 0.80) %>% 
    pull(ensembl_transcript_id)
  
# filter df_transcripts for these
  df_transcript_raw_counts_filtered <- df_transcript_raw_counts %>% 
    filter(ensembl_transcript_id %in% keep) # 70k

  

# normalize transcript counts using VST from DESEq2 ---------------------------------------------

  m_counts <- df_transcript_raw_counts_filtered %>% 
    as.data.frame %>% 
    column_to_rownames("ensembl_transcript_id") %>% 
    as.matrix + 1 # add a pseudo-count of 1 to avoid zeros in count mat
  
  m_vsd <- varianceStabilizingTransformation(m_counts, blind = TRUE, fitType = "parametric")



# map transcripts to DegTx using version id from biomart ------------------

  
# PULL TRANSCRIPT INFO FROM BIOMART  
  ensembl <- useEnsembl(biomart = "genes", 
                        dataset = "hsapiens_gene_ensembl")

  df_bm <- getBM(attributes = c("ensembl_transcript_id", "ensembl_transcript_id_version", "external_transcript_name", 
                                "ensembl_gene_id", "external_gene_name", 
                                "chromosome_name", "strand", 
                                "transcript_start", "transcript_end", "transcript_length", "transcript_biotype", "percentage_gene_gc_content"),
                 filters = "ensembl_transcript_id",
                 values = rownames(m_vsd),
                 mart = ensembl) %>% 
    as_tibble()
  
# GET DEGRADATION TRANSCRIPTS
  degradation_transcripts <- select_transcripts('cell_component')
  df_bm %>% filter(ensembl_transcript_id_version %in% degradation_transcripts) # 2,330
  
# MAP TRANSCRIPT ID VERSION TO TRANSCRIPT IN COUNT MATRIX
  m_vsd_filt <- m_vsd %>% 
    as.data.frame %>% 
    rownames_to_column("ensembl_transcript_id") %>% 
    as_tibble %>% 
    left_join(df_bm[1:2]) %>% 
    dplyr::select(ensembl_transcript_id_version, everything(), -ensembl_transcript_id) %>% 
    na.omit %>% 
    as.data.frame %>% 
    column_to_rownames("ensembl_transcript_id_version") %>% 
    as.matrix
  
  
  
  
# create RSE object --------------------------------------------------------------------



# ROW RANGES OBJECT
  gr <- GRanges(
    
    seqnames = Rle(df_bm$chromosome_name),
    ranges = IRanges(start = df_bm$transcript_start,
                     end = df_bm$transcript_end,
                     #width = df_bm$transcript_length,
                     names = rownames(m_vsd_filt)),
    strand = Rle(strand(df_bm$strand)),
    
    transcript_name = df_bm$external_transcript_name,
    type = df_bm$transcript_biotype,
    ensembl_gene_id = df_bm$ensembl_gene_id,
    gene_name = df_bm$external_gene_name,
    GC = df_bm$percentage_gene_gc_content
    
  )

# CREATE RSE
  rse_tx <- SummarizedExperiment(assays = SimpleList(counts = m_vsd_filt),
                                 rowRanges = gr,
                                 colData = m_metadata)

  

# run qSVA ----------------------------------------------------------------

  
# IDENTIFY DEGRADATION MATRIX FOR OUR TRANSCRIPTS
  DegTx <- getDegTx(rse_tx, 
                    type = "cell_component",
                    assayname = "counts")
  
# GET PCS FROM THIS MATRIX  
  pcTx <- getPCs(rse_tx = DegTx, 
                 assayname = "counts")

# DESIGN A BASIC MODEL MATRIX TO MODEL THE NUMBER OF PCS NEEDED  
  mod <- model.matrix(~ multi_nomial_dx + age_death + gender + race,
                      data = colData(rse_tx)
  )
  
# IDENTIFY NUMBER OF PCS NEEDED  
  set.seed(123)
  k <- k_qsvs(rse_tx = DegTx, 
              mod = mod, 
              assayname = "counts")
  
  print(k)

# USE K TO SUBSET OUR PCS
  qsvs <- get_qsvs(qsvPCs = pcTx, k = k)


# PERCENT VARIANCE OF DEGRADATION EXPLAINED BY qSVs
df_qsv_var <- summary(pcTx) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()

df_qsv_var %>% 
  filter(pc <= 16) %>% 
  
  ggplot(aes(x = pc)) +
  geom_col(aes(y = proportion_of_variance), fill = "midnightblue") +
  geom_point(aes(y = cumulative_proportion)) +
  geom_line(aes(y = cumulative_proportion)) +
  labs(x = "qSV") +
  ggtitle("proportion of degradation transcript matrix variance explained by qSVs 1-16")
  
df_qsv_var %>% 
  write.csv(paste0(base_dir, "outputs/qSV_explained_variance.csv"),
            row.names = FALSE)



# gene data ---------------------------------------------------------------

  
# COUNTS
  df_raw_counts <- read.table("data/185samples_allGenes_auto-PARs_noSexChrs_noFilters_KoryAnalysis_21Kgenes.txt") %>% 
    as_tibble(rownames = "ensembl_gene_id") %>% 
    rename_at(2:ncol(.), ~substr(.x, start = 2, stop = nchar(.x))) %>% 
    clean_names() %>% 
    dplyr::rename_at(2:ncol(.), ~str_remove(.x, "x"))

  

# remove low count genes --------------------------------------------------

  
  # 10 counts in at least 80% of samples
  keep <- df_raw_counts %>% 
    mutate_if(is.numeric, ~ifelse(.x > 10, 1, 0)) %>% 
    mutate(thresh = rowSums(dplyr::select(., -ensembl_gene_id))/(ncol(.) - 1), .before = 2) %>% 
    filter(thresh > 0.80) %>% 
    pull(ensembl_gene_id)
  
  # filter df_transcripts for these
  df_raw_counts_filtered <- df_raw_counts %>% 
    filter(ensembl_gene_id %in% keep) # 19k

  

# normalize counts using VST ----------------------------------------------

  
  m_counts <- df_raw_counts_filtered %>% 
  as.data.frame %>% 
  column_to_rownames("ensembl_gene_id") %>% 
  as.matrix + 1 # add a pseudo-count of 1 to avoid zeros in count mat

  m_vsd <- varianceStabilizingTransformation(m_counts, blind = TRUE, fitType = "parametric")
  
  df_vsd <- m_vsd %>% 
    as.data.frame %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    as_tibble

# COMBINE qSVS WITH COVARIATE DATA  
  df_covar_sv <- df_covariates_all_clean_numeric %>% 
    left_join(qsvs %>% 
                as.data.frame() %>% 
                rownames_to_column("sample"))
  

# PCA1: on normalized counts with no regression ---------------------------------------------------------------------
  
  
# PCA
  
  # PCA ON NORMALIZED COUNTS    
  
  pca <- prcomp(df_vsd %>%
                  as.data.frame %>% 
                  column_to_rownames("ensembl_gene_id"),
                center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    left_join(df_covar_sv)
  
  df_var <- summary(pca) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()
  
  # PLOT
  p_pca <- df_pca %>% 
    ggplot(aes(x = pc1, y = pc2)) +
    geom_point(color = "black", size = 2.75) +
    geom_point(aes(color = rin_acsg, text = sample), size = 2) +
    
    # label samples
    # geom_text_repel(aes(label = sample)) +
    
    # color gradient for continuous color variables
    scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
    
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
    ggtitle("PCA1: normalize counts")
  
  # plot eigenvalues
  p_var <- df_var %>% 
    ggplot(aes(x = pc, y = proportion_of_variance)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 0.01), lty = 2, color = "black")
  
  p_pca + p_var
  
  
  # CORRELATE PCs WITH COVARIATES
  df_pc_covar_cor <- df_pca %>% 
    pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
    pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
    group_by(pc, covariate) %>% 
    nest() %>% 
    
    mutate(pearsons_r = map(.x = data,
                            .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate
    ),
    p_val = map(.x = data,
                .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value
    )) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    ungroup %>% 
    mutate(q_val = p.adjust(p_val, method = "BH"))
  
  df_pc_covar_cor %>% 
    left_join(df_var %>% mutate(pc = paste0("pc", pc))) %>% 
    filter(proportion_of_variance > 0.01) %>% 
    mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
    mutate(labels = case_when(
      
      q_val < 0.05 & q_val > 0.01 ~ "*",
      q_val < 0.01 & q_val > 0.001 ~ "**",
      q_val < 0.001 ~ "***"
      
    )) %>% 
    
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = -log10(q_val))) +
    geom_text(aes(label = labels), color = "white", size = 8) +
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
    labs(x = "PC (> 1% explained var)") +
    theme_classic() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "bottom") +
    ggtitle("PCA1 correlations with covariates")
  
  
  
  
  
# PCA2: vsd ~ qSV1 + qSV4 + qSV6 + qSV13 + qSV15 ---------------------------------------------------------------------

# REGRESS ONLY qSV1 + qSV4 + qSV6 + qSV13 + qSV15
  
  df_vsd_regress <- df_vsd %>% 
    pivot_longer(2:ncol(.), names_to = "sample", values_to = "vsd") %>% 
    group_by(ensembl_gene_id) %>% 
    left_join(df_covar_sv) %>% 
    nest() %>% 
    
    mutate(resids = map(.x = data,
                        .f = ~ lm(vsd ~ qSV1 + qSV4 + qSV6 + qSV13 + qSV15, 
                                  data = .x %>% dplyr::select(vsd, contains("qSV"))) %>% 
                          residuals())) %>% 
    unnest(cols = c(data, resids))
  

# PCA2

  # PCA ON NORMALIZED & REGRESSED COUNTS    
  pca <- prcomp(df_vsd_regress %>%
                  pivot_wider(id_cols = ensembl_gene_id,
                              names_from = sample, values_from = resids) %>% 
                  as.data.frame %>% 
                  column_to_rownames("ensembl_gene_id"),
                center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    left_join(df_covar_sv)
    
  df_var <- summary(pca) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()
  
  # PLOT
  p_pca <- df_pca %>% 
    ggplot(aes(x = pc1, y = pc2)) +
    geom_point(color = "black", size = 2.75) +
    geom_point(aes(color = rin_acsg, text = sample), size = 2) +
    
    # label samples
    # geom_text_repel(aes(label = sample)) +
    
    # color gradient for continuous color variables
    scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
    
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
    ggtitle("PCA2: vst ~ qSV1 + qSV4 + qSV6 + qSV13 + qSV15")
  
  # plot eigenvalues
  p_var <- df_var %>% 
    ggplot(aes(x = pc, y = proportion_of_variance)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 0.01), lty = 2, color = "black")
  
  p_pca + p_var
  
  
# CORRELATE PCs WITH COVARIATES
  df_pc_covar_cor <- df_pca %>% 
    pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
    pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
    group_by(pc, covariate) %>% 
    nest() %>% 
    
    mutate(pearsons_r = map(.x = data,
                            .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate
    ),
    p_val = map(.x = data,
                .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value
    )) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    ungroup %>% 
    mutate(q_val = p.adjust(p_val, method = "BH"))
  
  df_pc_covar_cor %>% 
    left_join(df_var %>% mutate(pc = paste0("pc", pc))) %>% 
    filter(proportion_of_variance > 0.01) %>% 
    mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
    mutate(labels = case_when(
      
      q_val < 0.05 & q_val > 0.01 ~ "*",
      q_val < 0.01 & q_val > 0.001 ~ "**",
      q_val < 0.001 ~ "***"
      
    )) %>% 
    
    mutate(pc = str_remove(pc, "pc") %>% as.numeric %>% as.factor) %>% 
    
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = -log10(q_val))) +
    geom_text(aes(label = labels), color = "black", size = 5) +
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
    labs(x = "PC (> 1% explained var)") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "bottom") +
    ggtitle("correlate PCs with covariates")
  

  
    
# PCA3: vsd ~ qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 + qSV11 + qSV12 + qSV13 + qSV15 ---------------------------------------------------------------------
  
  
  # REGRESS qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 + qSV11 + qSV12 + qSV13 + qSV15
  
  df_vsd_regress <- df_vsd %>% 
    pivot_longer(2:ncol(.), names_to = "sample", values_to = "vsd") %>% 
    group_by(ensembl_gene_id) %>% 
    left_join(df_covar_sv) %>% 
    nest() %>% 
    
    mutate(resids = map(.x = data,
                        .f = ~ lm(vsd ~ qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 + qSV11 + qSV12 + qSV13 + qSV15, 
                                  data = .x %>% dplyr::select(vsd, contains("qSV"))) %>% 
                          residuals())) %>% 
    unnest(cols = c(data, resids))
  
  
# PCA3
  
  # PCA ON NORMALIZED & REGRESSED COUNTS    
  pca <- prcomp(df_vsd_regress %>%
                  pivot_wider(id_cols = ensembl_gene_id,
                              names_from = sample, values_from = resids) %>% 
                  as.data.frame %>% 
                  column_to_rownames("ensembl_gene_id"),
                center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    left_join(df_covar_sv)
  
  df_var <- summary(pca) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()
  
# PLOT
  p_pca <- df_pca %>% 
    ggplot(aes(x = pc1, y = pc2)) +
    geom_point(color = "black", size = 2.75) +
    geom_point(aes(color = rin_acsg, text = sample), size = 2) +
    
    # label samples
    # geom_text_repel(aes(label = sample)) +
    
    # color gradient for continuous color variables
    scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
    
    # add origin lines
    # geom_vline(xintercept = 0, color = "black", lty = 2) +
    # geom_hline(yintercept = 0, color = "black", lty = 2) +
    
    # give % variance in axis labels
    xlab(paste0("PC1: ", 
                df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance")) +
    ylab(paste0("PC2: ", 
                df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance"))
  
  # plot eigenvalues
  p_var <- df_var %>% 
    ggplot(aes(x = pc, y = proportion_of_variance)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 0.01), lty = 2, color = "black")
  
  p_pca + p_var + plot_annotation(title = "PCA3: vst ~ qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 + qSV11 + qSV12 + qSV13 + qSV15")
  
  
# CORRELATE PCs WITH COVARIATES
  df_pc_covar_cor <- df_pca %>% 
    pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
    pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
    group_by(pc, covariate) %>% 
    nest() %>% 
    
    mutate(pearsons_r = map(.x = data,
                            .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate
    ),
    p_val = map(.x = data,
                .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value
    )) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    ungroup %>% 
    mutate(q_val = p.adjust(p_val, method = "BH"))
  
  df_pc_covar_cor %>% 
    left_join(df_var %>% mutate(pc = paste0("pc", pc))) %>% 
    filter(proportion_of_variance > 0.01) %>% 
    mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
    mutate(labels = case_when(
      
      q_val < 0.05 & q_val > 0.01 ~ "*",
      q_val < 0.01 & q_val > 0.001 ~ "**",
      q_val < 0.001 ~ "***"
      
    )) %>% 
    
    mutate(pc = str_remove(pc, "pc") %>% as.numeric %>% as.factor) %>% 
    
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = -log10(q_val))) +
    geom_text(aes(label = labels), color = "black", size = 5) +
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
    labs(x = "PC (> 1% explained var)") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "bottom") +
    ggtitle("correlate PCs with covariates")
  
  

    
# PCA4: vsd ~ qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 + qSV9 + qSV10 + qSV12 + qSV13 + qSV14 + qSV15 ---------------------------------------------------------------------
  
  
  # REGRESS qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV9 + qSV8 + qSV10 + qSV12 + qSV13 + qSV14 + qSV15
  
  df_vsd_regress <- df_vsd %>% 
    pivot_longer(2:ncol(.), names_to = "sample", values_to = "vsd") %>% 
    group_by(ensembl_gene_id) %>% 
    left_join(df_covar_sv) %>% 
    nest() %>% 
    
    mutate(resids = map(.x = data,
                        .f = ~ lm(vsd ~ qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV9 + qSV8 + qSV10 + qSV12 + qSV13 + qSV14 + qSV15, 
                                  data = .x %>% dplyr::select(vsd, contains("qSV"))) %>% 
                          residuals())) %>% 
    unnest(cols = c(data, resids))
  
  
# PCA4
  
  # PCA ON NORMALIZED & REGRESSED COUNTS    
  pca <- prcomp(df_vsd_regress %>%
                  pivot_wider(id_cols = ensembl_gene_id,
                              names_from = sample, values_from = resids) %>% 
                  as.data.frame %>% 
                  column_to_rownames("ensembl_gene_id"),
                center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    left_join(df_covar_sv)
  
  df_var <- summary(pca) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()
  
  # PLOT
  p_pca <- df_pca %>% 
    ggplot(aes(x = pc1, y = pc2)) +
    geom_point(color = "black", size = 2.75) +
    geom_point(aes(color = rin_acsg, text = sample), size = 2) +
    
    # label samples
    # geom_text_repel(aes(label = sample)) +
    
    # color gradient for continuous color variables
    scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
    
    # add origin lines
    # geom_vline(xintercept = 0, color = "black", lty = 2) +
    # geom_hline(yintercept = 0, color = "black", lty = 2) +
    
    # give % variance in axis labels
    xlab(paste0("PC1: ", 
                df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance")) +
    ylab(paste0("PC2: ", 
                df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance"))
  
  # plot eigenvalues
  p_var <- df_var %>% 
    ggplot(aes(x = pc, y = proportion_of_variance)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 0.01), lty = 2, color = "black")
  
  p_pca + p_var + plot_annotation(title = "PCA4: vst ~ qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV9 + qSV8 + qSV10 + qSV12 + qSV13 + qSV14 + qSV15")
  
  
  # CORRELATE PCs WITH COVARIATES
  df_pc_covar_cor <- df_pca %>% 
    pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
    pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
    group_by(pc, covariate) %>% 
    nest() %>% 
    
    mutate(pearsons_r = map(.x = data,
                            .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate
    ),
    p_val = map(.x = data,
                .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value
    )) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    ungroup %>% 
    mutate(q_val = p.adjust(p_val, method = "BH"))
  
  df_pc_covar_cor %>% 
    left_join(df_var %>% mutate(pc = paste0("pc", pc))) %>% 
    filter(proportion_of_variance > 0.01) %>% 
    mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
    mutate(labels = case_when(
      
      q_val < 0.05 & q_val > 0.01 ~ "*",
      q_val < 0.01 & q_val > 0.001 ~ "**",
      q_val < 0.001 ~ "***"
      
    )) %>% 
    
    mutate(pc = str_remove(pc, "pc") %>% as.numeric %>% as.factor) %>% 
    
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = -log10(q_val))) +
    geom_text(aes(label = labels), color = "black", size = 5) +
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
    labs(x = "PC (> 1% explained var)") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "bottom") +
    ggtitle("correlate PCs with covariates")
  

    
# PCA5: vsd ~ qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 + qSV9 + qSV10 + qSV11 + qSV12 + qSV13 + qSV14 + qSV15 ---------------------------------------------------------------------
  
  
  # REGRESS qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 + qSV9 + qSV10 + qSV11 + qSV12 + qSV13 + qSV14 + qSV15
  
  df_vsd_regress <- df_vsd %>% 
    pivot_longer(2:ncol(.), names_to = "sample", values_to = "vsd") %>% 
    group_by(ensembl_gene_id) %>% 
    left_join(df_covar_sv) %>% 
    nest() %>% 
    
    mutate(resids = map(.x = data,
                        .f = ~ lm(vsd ~ qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 + qSV9 + qSV10 + qSV11 + qSV12 + qSV13 + qSV14 + qSV15, 
                                  data = .x %>% dplyr::select(vsd, contains("qSV"))) %>% 
                          residuals())) %>% 
    unnest(cols = c(data, resids))
  
  
# PCA5
  
  # PCA ON NORMALIZED & REGRESSED COUNTS    
  pca <- prcomp(df_vsd_regress %>%
                  pivot_wider(id_cols = ensembl_gene_id,
                              names_from = sample, values_from = resids) %>% 
                  as.data.frame %>% 
                  column_to_rownames("ensembl_gene_id"),
                center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    left_join(df_covar_sv)
  
  df_var <- summary(pca) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()
  
  # PLOT
  p_pca <- df_pca %>% 
    ggplot(aes(x = pc1, y = pc2)) +
    geom_point(color = "black", size = 2.75) +
    geom_point(aes(color = rin_acsg, text = sample), size = 2) +
    
    # label samples
    # geom_text_repel(aes(label = sample)) +
    
    # color gradient for continuous color variables
    scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
    
    # add origin lines
    # geom_vline(xintercept = 0, color = "black", lty = 2) +
    # geom_hline(yintercept = 0, color = "black", lty = 2) +
    
    # give % variance in axis labels
    xlab(paste0("PC1: ", 
                df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance")) +
    ylab(paste0("PC2: ", 
                df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance"))
  
  # plot eigenvalues
  p_var <- df_var %>% 
    ggplot(aes(x = pc, y = proportion_of_variance)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 0.01), lty = 2, color = "black")
  
  p_pca + p_var + 
    plot_annotation(title = "PCA5: vst ~ qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 + qSV9 + qSV10 + qSV11 + qSV12 + qSV13 + qSV14 + qSV15")
  
  
# CORRELATE PCs WITH COVARIATES
  df_pc_covar_cor <- df_pca %>% 
    pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
    pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
    group_by(pc, covariate) %>% 
    nest() %>% 
    
    mutate(pearsons_r = map(.x = data,
                            .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate
    ),
    p_val = map(.x = data,
                .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value
    )) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    ungroup %>% 
    mutate(q_val = p.adjust(p_val, method = "BH"))
  
  df_pc_covar_cor %>% 
    left_join(df_var %>% mutate(pc = paste0("pc", pc))) %>% 
    filter(proportion_of_variance > 0.01) %>% 
    mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
    mutate(labels = case_when(
      
      q_val < 0.05 & q_val > 0.01 ~ "*",
      q_val < 0.01 & q_val > 0.001 ~ "**",
      q_val < 0.001 ~ "***"
      
    )) %>% 
    
    mutate(pc = str_remove(pc, "pc") %>% as.numeric %>% as.factor) %>% 
    
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = -log10(q_val))) +
    geom_text(aes(label = labels), color = "black", size = 5) +
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
    labs(x = "PC (> 1% explained var)") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "bottom") +
    ggtitle("correlate PCs with covariates")
  
  
  
# PCA6: vsd ~ qSV1-16 ---------------------------------------------------------------------
  
  
  # REGRESS qSV1-16 
  
  df_vsd_regress <- df_vsd %>% 
    pivot_longer(2:ncol(.), names_to = "sample", values_to = "vsd") %>% 
    group_by(ensembl_gene_id) %>% 
    left_join(df_covar_sv) %>% 
    nest() %>% 
    
    mutate(resids = map(.x = data,
                        .f = ~ lm(vsd ~ ., 
                                  data = .x %>% dplyr::select(vsd, contains("qSV"))) %>% 
                          residuals())) %>% 
    unnest(cols = c(data, resids))
  
  
# PCA6
  
  # PCA ON NORMALIZED & REGRESSED COUNTS    
  pca <- prcomp(df_vsd_regress %>%
                  pivot_wider(id_cols = ensembl_gene_id,
                              names_from = sample, values_from = resids) %>% 
                  as.data.frame %>% 
                  column_to_rownames("ensembl_gene_id"),
                center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    left_join(df_covar_sv)
  
  df_var <- summary(pca) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()
  
  # PLOT
  p_pca <- df_pca %>% 
    ggplot(aes(x = pc1, y = pc2)) +
    geom_point(color = "black", size = 2.75) +
    geom_point(aes(color = rin_acsg, text = sample), size = 2) +
    
    # label samples
    # geom_text_repel(aes(label = sample)) +
    
    # color gradient for continuous color variables
    scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
    
    # add origin lines
    # geom_vline(xintercept = 0, color = "black", lty = 2) +
    # geom_hline(yintercept = 0, color = "black", lty = 2) +
    
    # give % variance in axis labels
    xlab(paste0("PC1: ", 
                df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance")) +
    ylab(paste0("PC2: ", 
                df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance"))
  
  # plot eigenvalues
  p_var <- df_var %>% 
    ggplot(aes(x = pc, y = proportion_of_variance)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 0.01), lty = 2, color = "black")
  
  p_pca + p_var + 
    plot_annotation(title = "PCA6: vst ~ qSV1-16")
  
  
# CORRELATE PCs WITH COVARIATES
  df_pc_covar_cor <- df_pca %>% 
    pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
    pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
    group_by(pc, covariate) %>% 
    nest() %>% 
    
    mutate(pearsons_r = map(.x = data,
                            .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate
    ),
    p_val = map(.x = data,
                .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value
    )) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    ungroup %>% 
    mutate(q_val = p.adjust(p_val, method = "BH"))
  
  df_pc_covar_cor %>% 
    left_join(df_var %>% mutate(pc = paste0("pc", pc))) %>% 
    filter(proportion_of_variance > 0.01) %>% 
    mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
    mutate(labels = case_when(
      
      q_val < 0.05 & q_val > 0.01 ~ "*",
      q_val < 0.01 & q_val > 0.001 ~ "**",
      q_val < 0.001 ~ "***"
      
    )) %>% 
    
    mutate(pc = str_remove(pc, "pc") %>% as.numeric %>% as.factor) %>% 
    
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = -log10(q_val))) +
    geom_text(aes(label = labels), color = "black", size = 5) +
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
    labs(x = "PC (> 1% explained var)") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "bottom") +
    ggtitle("correlate PCs with covariates")
  

  
# PCA7: vsd ~ qSV1-16 + mapped_percent---------------------------------------------------------------------
  
  
  # REGRESS qSV1-16 + mapped_percent
  
  df_vsd_regress <- df_vsd %>% 
    pivot_longer(2:ncol(.), names_to = "sample", values_to = "vsd") %>% 
    group_by(ensembl_gene_id) %>% 
    left_join(df_covar_sv) %>% 
    nest() %>% 
    
    mutate(resids = map(.x = data,
                        .f = ~ lm(vsd ~ ., 
                                  data = .x %>% dplyr::select(vsd, contains("qSV"), mapped_percent)) %>% 
                          residuals())) %>% 
    unnest(cols = c(data, resids))
  
  
# PCA7
  
  # PCA ON NORMALIZED & REGRESSED COUNTS    
  pca <- prcomp(df_vsd_regress %>%
                  pivot_wider(id_cols = ensembl_gene_id,
                              names_from = sample, values_from = resids) %>% 
                  as.data.frame %>% 
                  column_to_rownames("ensembl_gene_id"),
                center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    left_join(df_covar_sv)
  
  df_var <- summary(pca) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()
  
  # PLOT
  p_pca <- df_pca %>% 
    ggplot(aes(x = pc1, y = pc2)) +
    geom_point(color = "black", size = 2.75) +
    geom_point(aes(color = rin_acsg, text = sample), size = 2) +
    
    # label samples
    # geom_text_repel(aes(label = sample)) +
    
    # color gradient for continuous color variables
    scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
    
    # add origin lines
    # geom_vline(xintercept = 0, color = "black", lty = 2) +
    # geom_hline(yintercept = 0, color = "black", lty = 2) +
    
    # give % variance in axis labels
    xlab(paste0("PC1: ", 
                df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance")) +
    ylab(paste0("PC2: ", 
                df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance"))
  
  # plot eigenvalues
  p_var <- df_var %>% 
    ggplot(aes(x = pc, y = proportion_of_variance)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 0.01), lty = 2, color = "black")
  
  p_pca + p_var + 
    plot_annotation(title = "PCA7: vst ~ qSV1-16 + mapped_percent")
  
  
# CORRELATE PCs WITH COVARIATES
  df_pc_covar_cor <- df_pca %>% 
    pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
    pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
    group_by(pc, covariate) %>% 
    nest() %>% 
    
    mutate(pearsons_r = map(.x = data,
                            .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate
    ),
    p_val = map(.x = data,
                .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value
    )) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    ungroup %>% 
    mutate(q_val = p.adjust(p_val, method = "BH"))
  
  df_pc_covar_cor %>% 
    left_join(df_var %>% mutate(pc = paste0("pc", pc))) %>% 
    filter(proportion_of_variance > 0.01) %>% 
    mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
    mutate(labels = case_when(
      
      q_val < 0.05 & q_val > 0.01 ~ "*",
      q_val < 0.01 & q_val > 0.001 ~ "**",
      q_val < 0.001 ~ "***"
      
    )) %>% 
    
    mutate(pc = str_remove(pc, "pc") %>% as.numeric %>% as.factor) %>% 
    
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = -log10(q_val))) +
    geom_text(aes(label = labels), color = "black", size = 5) +
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
    labs(x = "PC (> 1% explained var)") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "bottom") +
    ggtitle("correlate PCs with covariates")
  

    
# PCA8: vsd ~ qSV1-16 + mapped_percent + rn_aextraction_batch ---------------------------------------------------------------------
  
  
  # REGRESS qSV1-16 + mapped_percent + rn_aextraction_batch
  
  df_vsd_regress <- df_vsd %>% 
    pivot_longer(2:ncol(.), names_to = "sample", values_to = "vsd") %>% 
    group_by(ensembl_gene_id) %>% 
    left_join(df_covar_sv) %>% 
    nest() %>% 
    
    mutate(resids = map(.x = data,
                        .f = ~ lm(vsd ~ ., 
                                  data = .x %>% dplyr::select(vsd, contains("qSV"), mapped_percent, rn_aextraction_batch)) %>% 
                          residuals())) %>% 
    unnest(cols = c(data, resids))
  
  
# PCA8
  
  # PCA ON NORMALIZED & REGRESSED COUNTS    
  pca <- prcomp(df_vsd_regress %>%
                  pivot_wider(id_cols = ensembl_gene_id,
                              names_from = sample, values_from = resids) %>% 
                  as.data.frame %>% 
                  column_to_rownames("ensembl_gene_id"),
                center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    left_join(df_covar_sv)
  
  df_var <- summary(pca) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()
  
  # PLOT
  p_pca <- df_pca %>% 
    ggplot(aes(x = pc1, y = pc2)) +
    geom_point(color = "black", size = 2.75) +
    geom_point(aes(color = rin_acsg, text = sample), size = 2) +
    
    # label samples
    # geom_text_repel(aes(label = sample)) +
    
    # color gradient for continuous color variables
    scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
    
    # add origin lines
    # geom_vline(xintercept = 0, color = "black", lty = 2) +
    # geom_hline(yintercept = 0, color = "black", lty = 2) +
    
    # give % variance in axis labels
    xlab(paste0("PC1: ", 
                df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance")) +
    ylab(paste0("PC2: ", 
                df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance"))
  
  # plot eigenvalues
  p_var <- df_var %>% 
    ggplot(aes(x = pc, y = proportion_of_variance)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 0.01), lty = 2, color = "black")
  
  p_pca + p_var + 
    plot_annotation(title = "PCA8: vst ~ qSV1-16 + mapped_percent + RNA_extraction_batch")
  
  
# CORRELATE PCs WITH COVARIATES
  df_pc_covar_cor <- df_pca %>% 
    pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
    pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
    group_by(pc, covariate) %>% 
    nest() %>% 
    
    mutate(pearsons_r = map(.x = data,
                            .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate
    ),
    p_val = map(.x = data,
                .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value
    )) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    ungroup %>% 
    mutate(q_val = p.adjust(p_val, method = "BH"))
  
  df_pc_covar_cor %>% 
    left_join(df_var %>% mutate(pc = paste0("pc", pc))) %>% 
    filter(proportion_of_variance > 0.01) %>% 
    mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
    mutate(labels = case_when(
      
      q_val < 0.05 & q_val > 0.01 ~ "*",
      q_val < 0.01 & q_val > 0.001 ~ "**",
      q_val < 0.001 ~ "***"
      
    )) %>% 
    
    mutate(pc = str_remove(pc, "pc") %>% as.numeric %>% as.factor) %>% 
    
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = -log10(q_val))) +
    geom_text(aes(label = labels), color = "black", size = 5) +
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
    labs(x = "PC (> 1% explained var)") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "bottom") +
    ggtitle("correlate PCs with covariates")
  
  



# correlate qSVs with known covariates ------------------------------------



df_qSV_covar_cor <- df_covar_sv %>% 
  pivot_longer(contains("qSV"), names_to = "qSV", values_to = "qSV_val") %>% 
  pivot_longer(!contains(c("qSV", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
  
  group_by(qSV, covariate) %>% 
  nest() %>% 
  
  mutate(pearsons_r = map(.x = data,
                          .f = ~ cor.test(.x$qSV_val, .x$covariate_val)$estimate
  ),
  p_val = map(.x = data,
              .f = ~ cor.test(.x$qSV_val, .x$covariate_val)$p.value
  )) %>% 
  unnest(cols = c(pearsons_r, p_val)) %>% 
  ungroup %>% 
  mutate(q_val = p.adjust(p_val, method = "BH"))

# PLOT
df_qSV_covar_cor %>% 
  mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
  mutate(labels = case_when(
    
    q_val < 0.05 & q_val > 0.01 ~ "*",
    q_val < 0.01 & q_val > 0.001 ~ "**",
    q_val < 0.001 ~ "***"
    
  )) %>% 
  
  mutate(qSV = str_remove(qSV, "qSV") %>% as.numeric %>% as.factor) %>% 
  
  ggplot(aes(x = qSV, y = covariate)) +
  geom_tile(aes(fill = -log10(q_val))) +
  geom_text(aes(label = labels), color = "black", size = 5) +
  scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
  labs(x = "qSV") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = "bottom") +
  ggtitle("correlate qSVs with known covariates")



# PCA9: vsd ~ qSV1-16 + mapped_percent + rn_aextraction_batch - qSV4 ---------------------------------------------------------------------



# REGRESS qSV1-16 + mapped_percent + rn_aextraction_batch - qSV4

df_vsd_regress <- df_vsd %>% 
  pivot_longer(2:ncol(.), names_to = "sample", values_to = "vsd") %>% 
  group_by(ensembl_gene_id) %>% 
  left_join(df_covar_sv) %>% 
  nest() %>% 
  
  mutate(resids = map(.x = data,
                      .f = ~ lm(vsd ~ . - qSV4, 
                                data = .x %>% dplyr::select(vsd, contains("qSV"), mapped_percent, rn_aextraction_batch)) %>% 
                        residuals())) %>% 
  unnest(cols = c(data, resids))


# PCA9

  # PCA ON NORMALIZED & REGRESSED COUNTS    
  pca <- prcomp(df_vsd_regress %>%
                  pivot_wider(id_cols = ensembl_gene_id,
                              names_from = sample, values_from = resids) %>% 
                  as.data.frame %>% 
                  column_to_rownames("ensembl_gene_id"),
                center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    left_join(df_covar_sv)
  
  df_var <- summary(pca) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()

  # PLOT
  p_pca <- df_pca %>% 
    ggplot(aes(x = pc1, y = pc2)) +
    geom_point(color = "black", size = 2.75) +
    geom_point(aes(color = rin_acsg, text = sample), size = 2) +
    
    # label samples
    # geom_text_repel(aes(label = sample)) +
    
    # color gradient for continuous color variables
    scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
    
    # add origin lines
    # geom_vline(xintercept = 0, color = "black", lty = 2) +
    # geom_hline(yintercept = 0, color = "black", lty = 2) +
    
    # give % variance in axis labels
    xlab(paste0("PC1: ", 
                df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance")) +
    ylab(paste0("PC2: ", 
                df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance"))

  # plot eigenvalues
  p_var <- df_var %>% 
    ggplot(aes(x = pc, y = proportion_of_variance)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 0.01), lty = 2, color = "black")
  
  p_pca + p_var + 
    plot_annotation(title = "PCA9: vst ~ qSV1-16 + mapped_percent + RNA_extraction_batch - qSV4")


# CORRELATE PCs WITH COVARIATES
df_pc_covar_cor <- df_pca %>% 
  pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
  pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
  group_by(pc, covariate) %>% 
  nest() %>% 
  
  mutate(pearsons_r = map(.x = data,
                          .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate
  ),
  p_val = map(.x = data,
              .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value
  )) %>% 
  unnest(cols = c(pearsons_r, p_val)) %>% 
  ungroup %>% 
  mutate(q_val = p.adjust(p_val, method = "BH"))

df_pc_covar_cor %>% 
  left_join(df_var %>% mutate(pc = paste0("pc", pc))) %>% 
  filter(proportion_of_variance > 0.01) %>% 
  mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
  mutate(labels = case_when(
    
    q_val < 0.05 & q_val > 0.01 ~ "*",
    q_val < 0.01 & q_val > 0.001 ~ "**",
    q_val < 0.001 ~ "***"
    
  )) %>% 
  
  mutate(pc = str_remove(pc, "pc") %>% as.numeric %>% as.factor) %>% 
  
  ggplot(aes(x = pc, y = covariate)) +
  geom_tile(aes(fill = -log10(q_val))) +
  geom_text(aes(label = labels), color = "black", size = 5) +
  scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
  labs(x = "PC (> 1% explained var)") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = "bottom") +
  ggtitle("correlate PCs with covariates")



# correlate qSVs with original gene-level modules -------------------------

  
  
# ORIGINAL RESIDUALIZED AND FILTERED GENE COUNTS
  df_vsd_resids <- read.csv("outputs/20220613_185samples_14kgenes_vsd_RIN_mappedPerc_GC_race_resids.csv") %>%
    as_tibble() %>%
    dplyr::select(-X) %>%
    rename_at(2:ncol(.), ~str_remove(.x, "X"))
  
# TRANPOSE COUNT RESIDUALS
  df_vsd_resids_t <- df_vsd_resids %>%
    
    # transpose & convert to tibble
    dplyr::select(-1) %>%
    t %>%
    as.data.frame %>%
    as_tibble %>%
    
    # rename columns with gene id from original df
    rename_all(~df_vsd_resids %>% pull(ensembl_gene_id)) %>%
    
    # add sample id col using colnames of original df (to maintain order)
    dplyr::mutate(sample = df_vsd_resids %>% dplyr::select(-1) %>% colnames, .before = 1)
  
  m_vsd_resids_t <- df_vsd_resids_t %>% as.data.frame %>% column_to_rownames("sample")
  
# LOAD MODULE ASSIGNMENTS  
  load("objects/archived/case_control/20220613_185samples_14kgenes_vsd_RIN_mappedPerc_GC_race_resids_sft10_11_minSize_cutHeight_mods.RDS")
  df_modules_filt <- df_modules %>% 
    filter(mod_set == "sft10_minSize30_cutHeight0.15") %>% 
    dplyr::select(-min_size, -cut_height)
  
# PULL MERGED MODULE EIGENGENES
  
  merged_colors <- df_modules_filt$merged_mod_color
  
  # assign colors to numbers
  df_color_mod <- df_modules_filt %>%
    dplyr::count(merged_mod_color) %>%
    arrange(desc(merged_mod_color == "grey"), desc(n)) %>%
    dplyr::rename("color" = "merged_mod_color") %>%
    mutate(module = 0:(n() - 1) %>% factor(levels = 0:(n() - 1)), .before = 1)
  
# GET EIGENGENES
  
  MEs <- moduleEigengenes(m_vsd_resids_t, merged_colors)$eigengenes
  
  df_MEs <- MEs %>%
    as.data.frame %>%
    rename_all(~str_remove(.x, "ME")) %>%
    rownames_to_column("sample") %>%
    as_tibble()
  
  df_MEs_long <- df_MEs %>%
    pivot_longer(2:ncol(.), names_to = "color", values_to = "kme") %>%
    left_join(df_color_mod, by = "color") %>%
    left_join(df_covariates_all_clean, by = "sample") %>%
    mutate(new_dx = factor(new_dx, levels = c("control", "bipolar", "mdd", "schizo")))
  
  df_sv_mod_cor <- df_MEs_long %>% 
    dplyr::select(sample, module, color, kme) %>% 
    left_join(df_covar_sv) %>% 
    dplyr::select(sample, module, color, kme, contains("SV")) %>% 
    pivot_longer(contains("SV"), names_to = "qSV_no", values_to = "value") %>% 
    
    group_by(module, color, qSV_no) %>% 
    nest() %>% 
    arrange(module) %>% 
    
    mutate(sv_module_r = map(.x = data,
                             .f = ~ cor.test(.x$kme, .x$value)$estimate),
           p_val = map(.x = data,
                       .f = ~ cor.test(.x$kme, .x$value)$p.value)) %>% 
    unnest(cols = c(sv_module_r, p_val)) %>% 
    mutate(q_val = p.adjust(p_val, method = "BH")) %>% 
    mutate(qSV_no = str_remove(qSV_no, "qSV") %>% as.numeric %>% factor(levels = 1:max(.)))
  
# PLOT MODULE-qSV ASSOCIATIONS    
  df_sv_mod_cor %>% 
    mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
    mutate(labels = case_when(
      
      q_val < 0.05 & q_val > 0.01 ~ "*",
      q_val < 0.01 & q_val > 0.001 ~ "**",
      q_val < 0.001 ~ "***"
      
    )) %>% 
    
    ggplot(aes(x = module, y = qSV_no)) +
    geom_tile(aes(fill = -log10(q_val))) +
    geom_text(aes(label = labels), color = "black", size = 8) +
    # geom_vline(aes(xintercept = 1.5), lty = 2, color = "black") +
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
    labs(x = "WGCNA module", y = "qSV") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "bottom") +
    ggtitle("correlate original WGCNA module eigengene (correct for RIN + mappedPerc + GCpercent + race) with qSV")
  
  
  
# PCA10: vsd ~ qSV1-7 + mapped_percent + rn_aextraction_batch - qSV4 ---------------------------------------------------------------------
  
  
  
  # REGRESS qSV1 + qSV2 + qSV3 + qSV5 + qSV6 + qSV7 + mapped_percent + rn_aextraction_batch
  
  df_vsd_regress <- df_vsd %>% 
    pivot_longer(2:ncol(.), names_to = "sample", values_to = "vsd") %>% 
    group_by(ensembl_gene_id) %>% 
    left_join(df_covar_sv) %>% 
    nest() %>% 
    
    mutate(resids = map(.x = data,
                        .f = ~ lm(vsd ~ qSV1 + qSV2 + qSV3 + qSV5 + qSV6 + qSV7 + mapped_percent + rn_aextraction_batch, 
                                  data = .x) %>% 
                          residuals())) %>% 
    unnest(cols = c(data, resids))
  
  
# PCA10
  
  # PCA ON NORMALIZED & REGRESSED COUNTS    
  pca <- prcomp(df_vsd_regress %>%
                  pivot_wider(id_cols = ensembl_gene_id,
                              names_from = sample, values_from = resids) %>% 
                  as.data.frame %>% 
                  column_to_rownames("ensembl_gene_id"),
                center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    left_join(df_covar_sv)
  
  df_var <- summary(pca) %>% 
    .$importance %>% 
    as.data.frame %>% 
    rownames_to_column("metric") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
    mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
    pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
    clean_names()
  
  # PLOT
  p_pca <- df_pca %>% 
    ggplot(aes(x = pc1, y = pc2)) +
    geom_point(color = "black", size = 2.75) +
    geom_point(aes(color = factor(gender), text = sample), size = 2) +
    # stat_conf_ellipse(aes(fill = factor(gender)), alpha = 0.5, geom = "polygon") +
    
    # label samples
    # geom_text_repel(aes(label = sample)) +
    
    # color gradient for continuous color variables
    #scale_color_viridis_b() +
    
    # add origin lines
    # geom_vline(xintercept = 0, color = "black", lty = 2) +
    # geom_hline(yintercept = 0, color = "black", lty = 2) +
    
    # give % variance in axis labels
    xlab(paste0("PC1: ", 
                df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance")) +
    ylab(paste0("PC2: ", 
                df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                "% of variance"))
  
  # plot eigenvalues
  p_var <- df_var %>% 
    ggplot(aes(x = pc, y = proportion_of_variance)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 0.01), lty = 2, color = "black")
  
  p_pca + p_var + 
    plot_annotation(title = "PCA: vst ~ qSV1-7 + mapped_percent + RNA_extraction_batch - qSV4")
  
  
  # CORRELATE PCs WITH COVARIATES
  df_pc_covar_cor <- df_pca %>% 
    pivot_longer(contains("pc"), names_to = "pc", values_to = "pc_val") %>% 
    pivot_longer(!contains(c("pc", "sample")), names_to = "covariate", values_to = "covariate_val") %>% 
    group_by(pc, covariate) %>% 
    nest() %>% 
    
    mutate(pearsons_r = map(.x = data,
                            .f = ~ cor.test(.x$pc_val, .x$covariate_val)$estimate
    ),
    p_val = map(.x = data,
                .f = ~ cor.test(.x$pc_val, .x$covariate_val)$p.value
    )) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    ungroup %>% 
    mutate(q_val = p.adjust(p_val, method = "BH"))
  
  df_pc_covar_cor %>% 
    left_join(df_var %>% mutate(pc = paste0("pc", pc))) %>% 
    filter(proportion_of_variance > 0.01) %>% 
    mutate(q_val = ifelse(p_val > 0.05, NA, q_val)) %>% 
    mutate(labels = case_when(
      
      q_val < 0.05 & q_val > 0.01 ~ "*",
      q_val < 0.01 & q_val > 0.001 ~ "**",
      q_val < 0.001 ~ "***"
      
    )) %>% 
    
    mutate(pc = str_remove(pc, "pc") %>% as.numeric %>% as.factor) %>% 
    
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = -log10(q_val))) +
    geom_text(aes(label = labels), color = "black", size = 5) +
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
    labs(x = "PC (> 1% explained var)") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "bottom") +
    ggtitle("correlate PCs with covariates")
  
  
base_dir <- "/Users/work/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/projects/WGCNA/acsg_wgcna/"
prefix <- "20230225_185samples_19kgenes_vst_qSVA123567_MP_RNAbatch_resids"
  
save(df_vsd_regress, file = paste0(base_dir, "objects/", prefix, ".RDS"))  
  
  
