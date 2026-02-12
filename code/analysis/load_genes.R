
#==============================================================================#
# Load gene-level data; map ensembl IDs to gene symbols
#==============================================================================#

## Note: This script can only be run after 00_preprocessing


base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "code/analysis/setup.R"))


### SET DIRECTORY AND LOAD DATA ###

prefix <- "08Mar2024_GENES_qSVAgeSexRaceGC"
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress


### MAP ENSEMBL IDs TO GENE SYMBOLS ###

ensembl_gene_ids <- df_vsd_regress %>% dplyr::select(-sample) %>% colnames
df_ensembl_to_symbol <- mapIds(
    EnsDb.Hsapiens.v79, 
    keys = ensembl_gene_ids,
    column = "SYMBOL", 
    keytype = "GENEID", 
    multiVals = "first"
) %>% 
    enframe(name = "ensembl_gene_id", value = "gene_symbol") %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol))
