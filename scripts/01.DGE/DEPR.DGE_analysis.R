
######################################################################

# Updated DGE analyis including all covariates used in WGCNA and GRCCA

######################################################################

## Set directories & load data
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/Fig2_DGEres/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_genes.R"))


# Convert count data to matrix for DESeq2; start by running only SCZ-control analyses ----------------------------

m_vsd_regress <- df_vsd_regress %>% 
    dplyr::filter(str_detect(sample, "schizo|control")) %>% # subset for only schizophrenia-control comparisons
    column_to_rownames("sample") %>% 
    t() + 10 # add pseudo-count of 10 so there are no negative values in the matrix
m_vsd_regress <- round(m_vsd_regress*10^5) # multiply by 10000 and round to convert to integers but maintain precision


# Generate covariates matrix ----------------------------------------------

## Combine drug MCs with known covariates
df_covariates_mca <- df_covariates %>% 
    left_join(
        df_ind_loadings %>% 
            dplyr::rename_at(vars(contains("dim")), ~ str_replace(.x, "dim", "MC")),
        by = join_by(sample)
    )

## Convert to matrix
n_mcs <- 8
m_covariates <- df_covariates_mca %>% 
    dplyr::filter(dx %in% c("SCZ", "Control")) %>% # SCZ vs controls only
    mutate(dx = factor(dx, levels = c("Control", "SCZ"))) %>% 
    dplyr::select(sample, dx, all_of(paste0("MC", 1:n_mcs)), brain_weight) %>% 
    column_to_rownames("sample")


# Run DESeq ---------------------------------------------------------------

## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = m_vsd_regress, 
                              colData = m_covariates, 
                              design = ~ dx + 
                                  MC1 + MC2 + MC3 + MC4 + MC5 + MC6 + MC7 + MC8 # + brain_weight; remove brain_weight only if regressing in CCA results (might be difficult to justify so removing for now)
                              )

## Run DESeq
de_res <- DESeq(dds)  

## Create tibble out of results
df_de_res <- results(de_res, name = "dx_SCZ_vs_Control") %>% 
    as.data.frame() %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    as_tibble() %>% 
    clean_names %>% 
    left_join(df_ensembl_to_symbol) %>% 
    dplyr::select(ensembl_gene_id, gene_symbol, everything()) %>% 
    arrange(-abs(log2fold_change)) %>% 
    mutate(padj = p.adjust(pvalue, method = "fdr")) %>% # overwrite the default DESeq2 padj because some of their values are NA
    mutate(base_mean = (base_mean/10^5) - 10) %>% 
    arrange(-abs(stat))


# Save for figures & downstream analysis -----------------------------------------

# Figures
save(df_de_res, file = paste0(analysis_objects_dir, "DE_results.RDS"))

# Downstream analyses
save(df_de_res, file = paste0(base_dir, "objects/DE_results.RDS"))

