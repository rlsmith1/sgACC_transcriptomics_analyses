
##############################################################################

# Run DESeq2 and Limma-voom without regressing technical covariates in advance
# (In order to fit model assumptions)

##############################################################################


# Setup -------------------------------------------------------------------


### SET DIRECTORIES ###

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/supplement/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")


### LOAD LIBRARIES ###

library(tidyverse)
library(DESeq2)
library(edgeR)


### LOAD DATA ###


## Gene raw counts
df_raw_counts <- read.table(paste0(base_dir, "data/185samples_allGenes_auto-PARs_noSexChrs_noFilters_KoryAnalysis_21Kgenes.txt")) %>% 
    as_tibble(rownames = "ensembl_gene_id") %>% 
    rename_at(2:ncol(.), ~substr(.x, start = 2, stop = nchar(.x))) %>% 
    clean_names() %>% 
    dplyr::rename_at(2:ncol(.), ~str_remove(.x, "x"))

## Cleaned covariates
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in 00.clean_covariates.R)

## Load QSVs (generated in 01.qSVA.R)
load(paste0(base_dir, "objects/06March2024_qSVs.Rdata"))

## Load MCA results
load(paste0(base_dir, "objects/drug_MCA_results.Rdata")) # df_ind_loadings, df_var_loadings



# Prepare covariate matrix ------------------------------------------------


## Combine qSVs and MCA with known data
df_covariates_all <- df_covariates %>% 
    left_join(df_qsvs) %>% 
    left_join(
        df_ind_loadings %>%
            rename_with(~ str_replace(.x, "dim", "MC"), contains("dim")),
        by = join_by(sample)
    )


## Define covariates that we want to include in downstream models
n_mcs <- 8
keep_qsvs <- paste0("qSV", c(seq(1, 7), 9)) # correlate with known covariates or var explained > 2% (see 02.gene_raw_count_preprocessing.R)
technical_covariates <- c("sex_at_birth", "race", "age_death", "gc_percent")


## Convert to matrix, keeping only covariates to include in the models
m_covariates <- df_covariates_all %>%
    dplyr::filter(dx %in% c("SCZ", "Control")) %>%
    
    # convert character variables to factors (for DESeq2)
    mutate(
        dx = factor(dx, levels = c("Control", "SCZ")),
        sex_at_birth = factor(sex_at_birth, levels = c("Female", "Male")),
        race = factor(race, levels = c("Black or African-American", "White"))
    ) %>%
    
    # select covariates to include
    dplyr::select(sample, dx, 
                  all_of( c(
                      paste0("MC", 1:n_mcs),
                      keep_qsvs,
                      technical_covariates
                  )
                  )
    ) %>%
    column_to_rownames("sample")



# Raw count filtering -----------------------------------------------------


## Filter for samples that are in gene data
df_raw_counts <- df_raw_counts %>% 
    dplyr::select(ensembl_gene_id, contains(df_covariates$sample))


## Remove low count genes

# At least 10 raw counts in at least 80% of samples
keep <- df_raw_counts %>% 
    mutate_if(is.numeric, ~ifelse(.x >= 10, 1, 0)) %>% 
    mutate(thresh = rowSums(dplyr::select(., -ensembl_gene_id))/(ncol(.) - 1), .before = 2) %>% 
    dplyr::filter(thresh >= 0.80) %>% 
    pull(ensembl_gene_id)
length(keep) # n = 18677

# filter df_raw_counts for genes that meet this threshold
df_raw_counts_filtered <- df_raw_counts %>% 
    dplyr::filter(ensembl_gene_id %in% keep)


## Match to samples of interest (SCZ and Control only)
df_raw_counts_filtered <- df_raw_counts_filtered %>% 
    dplyr::select(ensembl_gene_id, all_of(rownames(m_covariates)))


## Convert to matrix
m_counts <- df_raw_counts_filtered %>% 
    column_to_rownames("ensembl_gene_id")


## Check that sample ordering matches
all(colnames(m_counts) == rownames(m_covariates))


## Map ensembl IDs to gene symbols (save for later)
ensembl_gene_ids <- rownames(m_counts)
df_ensembl_to_symbol <- mapIds(
    EnsDb.Hsapiens.v79, 
    keys = ensembl_gene_ids,
    column = "SYMBOL", 
    keytype = "GENEID", 
    multiVals = "first"
) %>% 
    enframe(name = "ensembl_gene_id", value = "gene_symbol") %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol))



# DESeq2 ------------------------------------------------------------------


## Build a formula object (DESeq2 wants a formula, not a model matrix)
design_formula <- as.formula(paste0("~ ", paste0(colnames(m_covariates), collapse = " + ")))


## Construct DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = m_counts,
                              colData = m_covariates,
                              design = design_formula
)

## Run DESeq
dds <- DESeq(dds)

## Extract results
df_deseq2_res <- results(dds, name = "dx_SCZ_vs_Control") %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_gene_id") %>%
    as_tibble() %>%
    clean_names() %>%
    left_join(df_ensembl_to_symbol) %>%
    dplyr::select(ensembl_gene_id, gene_symbol, everything()) %>%
    arrange(-abs(log2fold_change)) %>%
    mutate(padj = p.adjust(pvalue, method = "fdr"))




# Limma-voom --------------------------------------------------------------


## Create DGEList
dge <- DGEList(counts = m_counts)
dge <- calcNormFactors(dge)

## Design matrix
design <- model.matrix(~ dx + ., data = m_covariates)

# ## Optional: filter genes again using edgeR’s method
# keep <- filterByExpr(dge, design)
# dge <- dge[keep, , keep.lib.sizes = FALSE]
# dge <- calcNormFactors(dge)  # Recalc after filtering

## Voom transformation
v <- voom(dge, design, plot = TRUE)

## Linear modeling
fit <- lmFit(v, design)
fit <- eBayes(fit)

## Extract results for dx
df_limma_res <- topTable(fit, coef = "dxSCZ", number = Inf) %>%
    rownames_to_column("ensembl_gene_id") %>%
    as_tibble() %>%
    clean_names() %>%
    left_join(df_ensembl_to_symbol) %>%
    dplyr::select(ensembl_gene_id, gene_symbol, everything()) %>%
    arrange(-abs(log_fc))




# Combine all differential expression results and save for plotting --------------------------------------------------------------


## Load Nirmala's DEG results
df_akula_res <- read_xlsx(path = paste0(base_dir, "data/Akula2021_geneLevel_21kgenes_baseMeanGe5_cqn-Rin-Race-GC_DESeq2-pvalsCorrected_betaPriorT-lfcShrink_092518.xlsx"),
                          sheet = "SczVsCtrls") %>% 
    clean_names %>% 
    dplyr::select(ensembl_gene_id, hgnc_symbol, base_mean, log2fold_change, lfc_se, stat, pvalue, padj) %>% 
    dplyr::rename("gene_symbol" = "hgnc_symbol")
    

## Combine DESeq2, limma, and Nirmala results
df_de_res_all <- df_deseq2_res %>% 
    dplyr::select(ensembl_gene_id, gene_symbol, log2fold_change, stat) %>%
    pivot_longer(3:4, names_to = "metric", values_to = "DESeq2") %>%
    mutate(metric = ifelse(str_detect(metric, "log"), "l2fc", "t_stat")) %>% 
    left_join(
        df_limma_res %>% 
            dplyr::select(ensembl_gene_id, gene_symbol, log_fc, t) %>% 
            pivot_longer(3:4, names_to = "metric", values_to = "Limma") %>% 
            mutate(metric = ifelse(str_detect(metric, "log"), "l2fc", "t_stat"))
    ) %>% 
    left_join(
        df_akula_res %>%
            dplyr::select(ensembl_gene_id, gene_symbol, log2fold_change, stat) %>% 
            pivot_longer(3:4, names_to = "metric", values_to = "Akula") %>% 
            mutate(metric = ifelse(str_detect(metric, "log"), "l2fc", "t_stat"))
    )


# Figures
save(
    df_deseq2_res, df_limma_res, df_akula_res, df_de_res_all, 
    file = paste0(analysis_objects_dir, "DE_results.Rdata")
)

# Downstream analyses
save(
    df_deseq2_res, df_limma_res, df_akula_res, df_de_res_all, 
    file = paste0(base_dir, "objects/DE_results.Rdata")
)




# Compare across methods ------------------------------------------------------------------



##
df_de_res_all %>% 
    ggplot(aes(x = Akula, y = DESeq2)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    geom_abline(lty = 2) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_cor() +
    facet_wrap(vars(metric), scales = "free")
    
##
df_de_res_all %>% 
    ggplot(aes(x = Limma, y = DESeq2)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    geom_abline(lty = 2) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_cor() +
    facet_wrap(vars(metric), scales = "free")

##
df_de_res_all %>% 
    ggplot(aes(x = Limma, y = Akula)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    geom_abline(lty = 2) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_cor() +
    facet_wrap(vars(metric), scales = "free")






