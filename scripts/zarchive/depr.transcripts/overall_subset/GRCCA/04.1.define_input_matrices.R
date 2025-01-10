
########################################################################################

# Define X.mat, Y.mat, and C.mat to include in GRCCA analysis

########################################################################################

### X.mat

# define gene order to align modules & genes
df_order <- df_modules_filt %>% 
    dplyr::select(ensembl_transcript_id, module) %>% 
    distinct() %>% 
    arrange(module, ensembl_transcript_id) #%>% 

## TURN OFF TO KEEP GRAY MODULE GENES
#filter(module != 0)

# X grouping vector    
x_group <- df_order %>% 
    mutate(module = as.numeric(as.character(module))) %>% 
    pull(module)

# X 
X.mat <- df_vsd_regress_filt %>% 
    dplyr::select(sample, all_of(df_order$ensembl_transcript_id)) %>% 
    column_to_rownames("sample")

sample_order <- rownames(X.mat)
transcripts <- colnames(X.mat)
length(transcripts) # n = 54302 (_WITHGRAY CVq1); # n = 24626 (_withGray CVq1)


### Y.MAT AND C.MAT

# define number of MCA dimensions to include
n_mca_dim <- 8

# Define variables to include in Y and C matrices
YC.mat <- df_covariates_numeric %>% 
    
    # Select covariates (MC loadings?)
    left_join(df_ind_loadings, by = join_by(sample)) %>%
    dplyr::select(sample, all_of(paste0("dim", seq(1, n_mca_dim)))
                  
    ) %>%
    mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
    
    # create dx columns
    mutate(control = ifelse(grepl("control", sample), 1, 0),
           bipolar = ifelse(grepl("bipolar", sample), 1, 0),
           mdd = ifelse(grepl("mdd", sample), 1, 0),
           schizo = ifelse(grepl("schizo", sample), 1, 0)
    ) %>% 
    dplyr::select(sample, bipolar, mdd, schizo, everything(), -control) %>% 
    mutate(sample = factor(sample, levels = sample_order)) %>% 
    arrange(sample) %>% 
    column_to_rownames("sample")

# Split into Y and C matrices
Y.mat <- YC.mat %>% dplyr::select(bipolar, mdd, schizo, all_of(paste0("dim", seq(1, n_mca_dim))), -dim2) # select suicide and drug MCs for Y.mat
head(Y.mat)
C.mat <- YC.mat %>% dplyr::select(-c(bipolar, mdd, schizo, all_of(paste0("dim", seq(1, n_mca_dim)))), dim2) # all other covariates go into C.mat
head(C.mat)

# check alignment of rows and columns
all(rownames(X.mat) == rownames(YC.mat))
all(colnames(X.mat) == c(df_order %>% pull(ensembl_transcript_id)))

### EXPORT MATRICES FOR CCA/PLS TOOLKIT

# define paths and create directory
type <- "OVERALLsubs/"
project_dir <- paste0(prefix,
                      "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "_", 
                      #"suicide", # is suicide included in Y.mat?
                      #"AllDrugs", # include all drugs in Y.mat
                      n_mca_dim, "MCA", # number of dimensions (if in Y.mat)
                      "_regress", # regressors included in C.mat
                      "Dim2", 
                      #n_mca_dim, "MCA", # number of dimensions (if in C.mat)
                      "_WITHGray/") # _WITHGRAY, _withGray
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)
dir.create(cca_dir, recursive = TRUE)
cca_data_dir <- paste0(cca_dir, "data/")
dir.create(cca_data_dir)

# write matrices
write.table(X.mat, file = paste0(cca_data_dir, "X.txt"))
write.table(Y.mat, file = paste0(cca_data_dir, "Y.txt"))
write.table(C.mat, file = paste0(cca_data_dir, "C.txt"))
write.table(x_group, file = paste0(cca_data_dir, "XGroup.txt"))

# GRCCA labels files
df_labels_x <- tibble(ensembl_transcript_id = colnames(X.mat)) %>% 
    left_join(df_modules_filt) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(Label = row_number(), .before = 1) %>% 
    dplyr::select(Label, module, ensembl_transcript_id, transcript_symbol) %>% 
    dplyr::rename_all(~c("Label", "Category", "Ensembl_id", "HGNC_id")) %>% 
    mutate(HGNC_id = ifelse(is.na(HGNC_id), Ensembl_id, HGNC_id))

df_labels_y <- tibble(Label = colnames(Y.mat))

write.csv(df_labels_x, paste0(cca_data_dir, "LabelsX.csv"), row.names = FALSE)
write.csv(df_labels_y, paste0(cca_data_dir, "LabelsY.csv"), row.names = FALSE)



