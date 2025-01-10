
######################################################################################

# Identify qSVs using transcript counts and qsvaR to regress in step 01
# (same qSVs for both gene and transcript level)

######################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(janitor)
library(qsvaR)
library(biomaRt)
library(RColorBrewer)
library(patchwork)
library(DESeq2)


# Load transcript & covariate data ---------------------------------------------------------


# TRANSCRIPT COUNTS
  df_transcript_raw_counts <- read.csv(paste0(base_dir, "data/transcript_count_matrix.csv")) %>% 
    as_tibble() %>% 
    dplyr::rename("ensembl_transcript_id" = "X") %>% 
    rename_all(~str_remove(.x, "_MergedBam_stringtieOutput")) %>% 
    clean_names() %>% 
    rename_all(~str_remove(.x, "x")) 

# COVARIATES
  load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_clean (generated in clean_covariates.R)
  m_metadata <- df_covariates_numeric %>% 
    as.data.frame %>% 
    column_to_rownames("sample") %>% 
    as.matrix

# FILTER COUNTS FOR ONLY SAMPLES THAT ARE IN GENE DATA
  df_transcript_raw_counts <- df_transcript_raw_counts %>% 
    dplyr::select(ensembl_transcript_id, contains(df_covariates_numeric$sample))
  
  
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

  
# Normalize transcript counts using VST from DESEq2 ---------------------------------------------

  m_counts <- df_transcript_raw_counts_filtered %>% 
    as.data.frame %>% 
    column_to_rownames("ensembl_transcript_id") %>% 
    as.matrix + 1 # add a pseudo-count of 1 to avoid zeros in count mat
  
  m_vsd <- varianceStabilizingTransformation(m_counts, blind = TRUE, fitType = "parametric")


# Map transcripts to degradation transcripts using version id from biomart ------------------
  
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
  df_bm %>% filter(ensembl_transcript_id_version %in% degradation_transcripts) # 2,292
  
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
  

  
# Create RSE object --------------------------------------------------------------------

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


# Run qSVA ----------------------------------------------------------------
  

# IDENTIFY DEGRADATION MATRIX FOR OUR TRANSCRIPTS
  DegTx <- getDegTx(rse_tx,
                    type = "cell_component",
                    assayname = "counts")

# GET PCS FROM THIS MATRIX (to plot variance explained by each qSV)
  set.seed(20240306)
  pcTx <- getPCs(rse_tx = DegTx,
                 assayname = "counts")
  
  df_qsv_var_expl <- summary(pcTx)$importance %>% 
    as.data.frame() %>% 
    rownames_to_column("feature") %>% 
    as_tibble %>% 
    pivot_longer(2:ncol(.), names_to = "PC", values_to = "value") %>% 
    mutate(qSV = str_remove(PC, "PC") %>% as.numeric %>% factor(levels = c(1:max(.)))) 

# DESIGN A BASIC MODEL MATRIX TO MODEL THE NUMBER OF PCS NEEDED  
  mod <- model.matrix(~ dx,
                      data = colData(rse_tx)
  )
  
# WRAPPER FUNCTION TO EXTRACT qSVs
  set.seed(20240306)
  qsvs <- qSVA(rse_tx = rse_tx, type = "cell_component", mod = mod, assayname = "counts")
  
# CONVERT TO TIBBLE AND SAVE
  df_qsvs <- qsvs %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble()
  

# Correlate qSVs with known covariates ------------------------------------

## Identify covariates and qSVs to correlate
  covariates <- df_covariates_numeric %>% dplyr::select(-sample) %>% colnames
  qsvs <- df_qsvs %>% dplyr::select(-sample) %>% colnames
  
## Run all combinations of correlations  
  df_covariate_qsv_cor <- expand_grid(
      covariate = covariates,
      qsv = qsvs
  ) %>% 
      mutate(
          pearsons_r = map2(
              .x = covariate,
              .y = qsv,
              .f = ~ cor.test(df_covariates_numeric[[.x]], df_qsvs[[.y]])$estimate
          ),
          p_value = map2(
              .x = covariate,
              .y = qsv,
              .f = ~ cor.test(df_covariates_numeric[[.x]], df_qsvs[[.y]])$p.value
          )
      ) %>% 
      unnest(cols = c(pearsons_r, p_value)) %>% 
      mutate(p_adj = p.adjust(p_value, method = "fdr")) %>% 
      mutate(qsv = str_remove(qsv, "qSV") %>% 
                 as.numeric %>% 
                 as.factor)
  
  

# Save objects for figures ------------------------------------------------

  save(df_qsvs, df_qsv_var_expl, df_covariate_qsv_cor,
       file = paste0(analysis_objects_dir, "06March2024_qSVs.Rdata")
  ) # df_qsvs, df_qsv_var_expl
  

