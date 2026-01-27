
#==============================================================================#
# Identify qSVs using transcript counts and qsvaR to regress in step 01
# (same qSVs for both gene and transcript level)
#==============================================================================#

## Run qSVA using qsvaR to account for transcript degradation in post-mortem
## gene expression data
## Will combine with known covariates for covariate selection


# Setup ---------------------------------------------------------

## Set directories
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")


## Libraries
library(tidyverse)
library(janitor)
library(qsvaR)
library(biomaRt)
library(RColorBrewer)
library(patchwork)
library(DESeq2)


## Load data

# Transcript count data
df_transcript_raw_counts <- read.csv(paste0(base_dir, "data/transcript_count_matrix.csv")) %>% 
    as_tibble() %>% 
    dplyr::rename("ensembl_transcript_id" = "X") %>% 
    rename_all(~str_remove(.x, "_MergedBam_stringtieOutput")) %>% 
    clean_names() %>% 
    rename_all(~str_remove(.x, "x")) 

# Formatted covariate data
load(paste0(base_dir, "objects/22Jan2025_covariates.Rdata")) # df_covariates, df_covariates_clean (generated in 00_clean-covariates.R)

# Convert to matrix
m_metadata <- df_covariates_numeric %>% 
    as.data.frame %>% 
    column_to_rownames("sample") %>% 
    as.matrix

# Filter count data for samples that are included in covariate data
df_transcript_raw_counts <- df_transcript_raw_counts %>% 
    dplyr::select(ensembl_transcript_id, contains(df_covariates_numeric$sample))


  
# Remove low count transcripts --------------------------------------------

## At least 10 raw counts in at least 80% of samples
keep <- df_transcript_raw_counts %>% 
    mutate_if(is.numeric, ~ifelse(.x >= 10, 1, 0)) %>% 
    mutate(thresh = rowSums(dplyr::select(., -ensembl_transcript_id))/(ncol(.) - 1), .before = 2) %>% 
    filter(thresh >= 0.80) %>% 
    pull(ensembl_transcript_id)
length(keep) # 72403
  
## Filter df_transcripts for transcripts that meet filtering criteria
df_transcript_raw_counts_filtered <- df_transcript_raw_counts %>% 
    filter(ensembl_transcript_id %in% keep)

  
# Normalize transcript counts using VST from DESEq2 ---------------------------------------------

m_counts <- df_transcript_raw_counts_filtered %>% 
    as.data.frame %>% 
    column_to_rownames("ensembl_transcript_id") %>% 
    as.matrix + 1 # add a pseudo-count of 1 to avoid zeros in count mat

m_vsd <- varianceStabilizingTransformation(m_counts, blind = TRUE, fitType = "parametric")


# Map transcripts to degradation transcripts using version id from biomart ------------------
  
# Use biomart to get transcript info 
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
  
## Get degradation transcripts from qsvaR
degradation_transcripts <- select_transcripts('cell_component')
df_bm %>% filter(ensembl_transcript_id_version %in% degradation_transcripts) # 2,292

## Map transcript ID version to transcript in normalized count matrix
m_vsd_filt <- m_vsd %>% 
    
    # convert to tibble
    as.data.frame %>% 
    rownames_to_column("ensembl_transcript_id") %>% 
    as_tibble %>% 
    
    # join with biomart info
    left_join(df_bm[1:2]) %>% 
    dplyr::select(ensembl_transcript_id_version, everything(), -ensembl_transcript_id) %>% 
    na.omit %>% 
    
    # convert back to matrix
    as.data.frame %>% 
    column_to_rownames("ensembl_transcript_id_version") %>% 
    as.matrix


  
# Create RSE object --------------------------------------------------------------------

## Row ranges object
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

## Create RSE (for qSVA input)
rse_tx <- SummarizedExperiment(assays = SimpleList(counts = m_vsd_filt),
                               rowRanges = gr,
                               colData = m_metadata)


# Run qSVA ----------------------------------------------------------------
  

## Identify degradation matrix for our transcripts
  DegTx <- getDegTx(rse_tx,
                    type = "cell_component",
                    assayname = "counts")

## Extract PCs from this matrix
set.seed(20240306)
pcTx <- getPCs(rse_tx = DegTx,
               assayname = "counts")

df_qsv_var_expl <- summary(pcTx)$importance %>% 
    as.data.frame() %>% 
    rownames_to_column("feature") %>% 
    as_tibble %>% 
    pivot_longer(2:ncol(.), names_to = "PC", values_to = "value") %>% 
    mutate(qSV = str_remove(PC, "PC") %>% as.numeric %>% factor(levels = c(1:max(.)))) 

## Design a basic model matrix to model the number of PCs (qSVs) needed
mod <- model.matrix(~ dx,
                    data = colData(rse_tx)
)

## Wrapper function to extract qSVs
set.seed(20240306)
qsvs <- qSVA(rse_tx = rse_tx, type = "cell_component", mod = mod, assayname = "counts")

## Convert to tibble
df_qsvs <- qsvs %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble()
  

# Correlate qSVs with known covariates ------------------------------------

## Identify covariates and qSVs to correlate
covariates <- df_covariates_numeric %>% dplyr::select(-sample) %>% colnames
qsvs <- df_qsvs %>% dplyr::select(-sample) %>% colnames
  
## Run all pairwise combinations of correlations  
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
  
  

# Save objects ------------------------------------------------

  save(df_qsvs, df_qsv_var_expl, df_covariate_qsv_cor,
       file = paste0(analysis_objects_dir, "06March2024_qSVs.Rdata")
  ) # df_qsvs, df_qsv_var_expl
  

