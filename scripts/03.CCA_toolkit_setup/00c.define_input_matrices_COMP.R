
######################################################################################

# Define X.mat, Y.mat, and C.mat to include in (G)RCCA analysis -- COMPOSITE CASE

######################################################################################


### Setup ###

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/load_WGCNA_res.R"))


### SELECT samples to include ###
dx_groups_to_analyze <- c("Control", "BD", "SCZ", "MDD") # select diagnostic groups of interest
samples_to_include <- df_covariates %>% 
    filter(dx %in% dx_groups_to_analyze) %>% 
    pull(sample) %>% 
    unique



### X.MAT ###

# define gene order to align modules & genes
df_order <- df_modules_filt %>% 
    dplyr::select(ensembl_gene_id, module) %>% 
    distinct() %>% 
    arrange(module, ensembl_gene_id) #%>% 
    
    # remove gray module!
    #filter(module != 0)

# X grouping vector    
x_group <- df_order %>% 
    mutate(module = as.numeric(as.character(module %>% str_remove("geneM")))) %>% 
    pull(module)

# X 
X.mat <- df_vsd_regress %>% 
    filter(sample %in% samples_to_include) %>% 
    dplyr::select(sample, all_of(df_order$ensembl_gene_id)) %>% 
    column_to_rownames("sample")
sample_order <- rownames(X.mat)

dim(X.mat) # n = 18677



##### COMPOSITE CASE: Combine SCZ, BD, and MDD into 1 diagnosis column (1 for psych diagnosis, 0 for no known diagnosis) #####

### Y.MAT AND C.MAT

# define number of MCA dimensions to include
n_mca_dim <- 8

# Define variables to include in Y and C matrices
YC.mat <- df_covariates_numeric %>% 
    
    # include only samples selected
    filter(sample %in% samples_to_include) %>% 
    
    # drugs to keep
    # dplyr::select(sample, age_death, brain_weight, opioids,
    #               sedative_hypnotic_anxiolitics, antidepressants, antipsychotics, mood_stabilizers, benzos) %>%
    # mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
    
    # Add MCA loadings
    left_join(df_ind_loadings, by = join_by(sample)) %>%
    
    # create dx column <------------------------------------ # DEFINE COMPOSITE CASE HERE #
    mutate(dx = ifelse(grepl("control", sample), 0, 1)) %>% 
    
    # select which variables to include in the Y & C (covariate and confound) matrices
    dplyr::select(sample, 
                  dx,
                  #all_of(dx_groups_to_analyze),
                  all_of(paste0("dim", seq(1, n_mca_dim))),
                  brain_weight
    ) %>%
    mutate(sample = factor(sample, levels = sample_order)) %>% 
    arrange(sample) %>% 
    column_to_rownames("sample")

# Split into Y and C matrices
Y.mat <- YC.mat %>% dplyr::select(-c(brain_weight))
head(Y.mat)
C.mat <- YC.mat %>% dplyr::select(brain_weight)
head(C.mat)

# check alignment of rows and columns
all(rownames(X.mat) == rownames(YC.mat))
all(colnames(X.mat) == c(df_order %>% pull(ensembl_gene_id)))



### EXPORT MATRICES FOR CCA/PLS TOOLKIT

# define paths and create directory
type <- "GENES/"
project_dir <- paste0("26Sept2025_COMPOSITE",
                      "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "/"
                      #"_seed", seed, "/" # includes gray module, does not regress brain weight
                      #paste0(dx_groups_to_analyze, collapse = ""), "_", # diagnostic groups included
                      #n_mca_dim, "MCA_regressBrainWeight_withGray/" # _WITHGRAY, _withGray
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
df_labels_x <- tibble(ensembl_gene_id = colnames(X.mat)) %>% 
    left_join(df_modules_filt) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(Label = row_number(), .before = 1) %>% 
    dplyr::select(Label, module, ensembl_gene_id, gene_symbol) %>% 
    dplyr::rename_all(~c("Label", "Category", "Ensembl_id", "HGNC_id")) %>% 
    mutate(HGNC_id = ifelse(is.na(HGNC_id), Ensembl_id, HGNC_id))

df_labels_y <- tibble(Label = colnames(Y.mat))

write.csv(df_labels_x, paste0(cca_data_dir, "LabelsX.csv"), row.names = FALSE)
write.csv(df_labels_y, paste0(cca_data_dir, "LabelsY.csv"), row.names = FALSE)

