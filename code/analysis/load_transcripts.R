
#==============================================================================#
# Load transcript-level data; map transcript IDs to genes
#==============================================================================#

## Note: This script can only be run after 00_preprocessing

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "code/analysis/setup.R"))


### SET DIRECTORY AND LOAD DATA ###

prefix <- "08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC_CVq1"
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress_filt


### MAP TRANSCRIPTS TO GENES ###

### Transcript ensembl ID to gene ID & symbol mappings ###
ensembl_transcript_ids <- df_vsd_regress_filt %>% dplyr::select(-sample) %>% colnames %>% unique

## Access biomaRt database to map transcript IDs to gene IDs
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
biomart_res <- getBM(attributes = c("ensembl_transcript_id", 
                                    "external_transcript_name",
                                    "ensembl_gene_id"),
                     filters = "ensembl_transcript_id", 
                     values = ensembl_transcript_ids,
                     mart = mart
) %>% 
    as_tibble() %>% 
    dplyr::rename("transcript_symbol" = "external_transcript_name")
rm(list = "mart")

## Identify gene symbol for genes that were not mapped
df_ensembl_to_symbol <- mapIds(
    EnsDb.Hsapiens.v79, 
    keys = biomart_res$ensembl_gene_id,
    column = "SYMBOL", 
    keytype = "GENEID", 
    multiVals = "first"
) %>% 
    enframe(name = "ensembl_gene_id", value = "gene_symbol") %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol))

## Combine
df_transcript_to_gene <- biomart_res %>% 
    left_join(df_ensembl_to_symbol, relationship = "many-to-many", by = join_by(ensembl_gene_id)) %>% 
    mutate(transcript_symbol = ifelse(transcript_symbol == "", NA_character_, transcript_symbol)) %>% 
    mutate(transcript_symbol = ifelse(is.na(transcript_symbol) & !str_detect(gene_symbol, "^ENSG"), gene_symbol, transcript_symbol))

