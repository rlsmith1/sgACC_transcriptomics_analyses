
################################################################################

# Define X.mat, Y.mat, and C.mat to include in GRCCA analysis

################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/load_WTCNA_res.R"))


### SELECT samples to include ###
dx_groups_to_analyze <- c("Control", "BD", "MDD", "SCZ") # select diagnostic groups of interest
samples_to_include <- df_covariates %>% 
    filter(dx %in% dx_groups_to_analyze) %>% 
    pull(sample) %>% 
    unique

### X.mat ###

# NEW: filter transcripts for risk genes derivatives and/or genes that were interesting in GRCCA
load(paste0(base_dir, "outputs/objects/Fig2_GRCCAres/GRCCA_results.Rdata")) # df_x_res
df_x_res_genes <- df_x_res
rm(list = "df_x_res")

grcca_genes_to_include <- df_x_res_genes %>% 
    filter(p_adj < 0.05) %>% # significantly associated with latent variable (not necessarily non-zero)
    #filter(abs(z_score) >= 1.96) %>% # genes with significant weights (not necessarily related to SCZ signal)
    pull(ensembl_gene_id) %>%
    unique
risk_genes_to_include <- unlist(benchmarking_lists) %>% unique # all psychiatric risk genes
genes_to_include <- unique(c(grcca_genes_to_include, risk_genes_to_include))
length(genes_to_include) # n = 5130 (FDR < 0.05 & risk genes); n = 2714 (FDR < 0.05 & |z| >= 2); n = 4491 (|z| >= 1.96)

transcripts_to_include <- df_transcript_to_gene %>% 
    filter(ensembl_gene_id %in% genes_to_include) %>% 
    pull(ensembl_transcript_id) %>% 
    unique()
length(transcripts_to_include) # n = 6705 (&); n = 1876 (broad & rare SCZ risk genes); 
# n = 12986 for GRCCA FDR < 0.05 and psychiatric risk genes
# n = 10344 (significant GRCCA genes FDR < 0.05); 
# n = 3610 (all psychiatric risk genes)
# 11040 for abs(z) >= 1.96
# 48213 for only genes included at gene-level

# define gene order to align modules & genes
df_order <- df_modules_filt %>% 
    filter(ensembl_transcript_id %in% transcripts_to_include) %>% #*filter here!!!*#
    dplyr::select(ensembl_transcript_id, module) %>% 
    distinct() %>% 
    arrange(module, ensembl_transcript_id)

# X grouping vector    
x_group <- df_order %>% 
    mutate(module = as.numeric(as.character(module %>% str_remove("transcriptM")))) %>% 
    pull(module)

# X 
X.mat <- df_vsd_regress_filt %>%
    filter(sample %in% samples_to_include) %>%
    dplyr::select(sample, all_of(df_order$ensembl_transcript_id)) %>% 
    column_to_rownames("sample")
dim(X.mat)
sample_order <- rownames(X.mat)
transcripts <- colnames(X.mat)
length(transcripts) # n = 54302 (for all transcripts)


### Y.MAT AND C.MAT

# define number of MCA dimensions to include
n_mca_dim <- 8

# Define variables to include in Y and C matrices
YC.mat <- df_covariates_numeric %>%
    
    # include only samples selected
    filter(sample %in% samples_to_include) %>% 
    
    # Select covariates (MC loadings?)
    left_join(df_ind_loadings, by = join_by(sample)) %>%
    #mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
    
    # create dx columns
    mutate(Control = ifelse(grepl("control", sample), 1, 0),
           BD = ifelse(grepl("bipolar", sample), 1, 0),
           MDD = ifelse(grepl("mdd", sample), 1, 0),
           SCZ = ifelse(grepl("schizo", sample), 1, 0)
    ) %>% 
    
    # select which variables to include in the Y & C (covariate and confound) matrices
    dplyr::select(sample, 
                  all_of(dx_groups_to_analyze),
                  all_of(paste0("dim", seq(1, n_mca_dim))),
                  #brain_weight,
                  -Control # remove control as matrix intercept
    ) %>%
    mutate(sample = factor(sample, levels = sample_order)) %>% 
    arrange(sample) %>% 
    column_to_rownames("sample")

# Split into Y and C matrices
head(YC.mat)
Y.mat <- YC.mat %>% dplyr::select(-c(dim2))
head(Y.mat)
C.mat <- YC.mat %>% dplyr::select(dim2)
head(C.mat)

# check alignment of rows and columns
all(rownames(X.mat) == rownames(YC.mat))
all(colnames(X.mat) == c(df_order %>% pull(ensembl_transcript_id)))


### EXPORT MATRICES FOR CCA/PLS TOOLKIT

# define paths and create directory
type <- "OVERALLsubs/"
project_dir <- paste0(prefix, # date & preprocessing info
                      "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "_", # WTCNA parameters
                      paste0(dx_groups_to_analyze, collapse = ""), "_", # diagnostic groups included
                      n_mca_dim, "MCA_regressDim2_", # covariates included/regressed
                      "broadRareSCZBDMDDASD_sigGRCCAfdr0_05/" #"sigGRCCAfdr0_05/" # define subset of transcripts included (based on gene-level results)
)
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)
dir.create(cca_dir, recursive = TRUE)
cca_data_dir <- paste0(cca_dir, "data/")
dir.create(cca_data_dir)

# write matrices
write.table(X.mat, file = paste0(cca_data_dir, "X.txt"))
write.table(Y.mat, file = paste0(cca_data_dir, "Y.txt"))
write.table(C.mat, file = paste0(cca_data_dir, "C.txt"))
write.table(x_group, file = paste0(cca_data_dir, "GroupsX.txt"))

# GRCCA labels files
df_labels_x <- tibble(ensembl_transcript_id = colnames(X.mat)) %>% 
    left_join(df_modules_filt) %>% 
    left_join(df_transcript_to_gene, by = join_by(ensembl_transcript_id)) %>%
    dplyr::select(module, ensembl_transcript_id, transcript_symbol) %>%
    distinct %>% 
    mutate(Label = row_number(), .before = 1) %>% 
    dplyr::rename_all(~c("Label", "Category", "Ensembl_id", "HGNC_id")) %>% 
    mutate(HGNC_id = ifelse(is.na(HGNC_id), Ensembl_id, HGNC_id))

df_labels_y <- tibble(Label = colnames(Y.mat))

write.csv(df_labels_x, paste0(cca_data_dir, "LabelsX.csv"), row.names = FALSE)
write.csv(df_labels_y, paste0(cca_data_dir, "LabelsY.csv"), row.names = FALSE)
