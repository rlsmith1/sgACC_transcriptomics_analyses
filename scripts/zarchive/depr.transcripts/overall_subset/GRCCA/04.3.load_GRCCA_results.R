
########################################################################################

# Load GRCCA results to analyze

########################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/full_analysis_scripts/transcripts/overall_subset/GRCCA/setup.R"))

# Load toolkit output -----------------------------------------------------

# READ IN GRCCA RES (GENERATED ON CLUSTER)
n_mca_dim <- 8
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
analysis_dir <- "grcca_permutation_VARx0.1_1_L2xNfeat_groupL2x0.1_withCmat_allEffects" #_noCmat 
model_1 <- read_mat(paste0(cca_dir, "framework/", analysis_dir, "/res/level1/model_1.mat"))
boot_1 <- read_mat(paste0(cca_dir, "framework/", analysis_dir, "/boot/level1/allboot_1.mat"))


# READ RESULTS TABLE TO SELECT SPLIT
df_results <- read_table(paste0(cca_dir, "framework/", analysis_dir, "/res/level1/results_table.txt"))
df_results
best_split <- df_results %>% arrange(pval, -correl) %>% pull(set) %>% .[[1]]


# Load X.mat & Y.mat ------------------------------------------------------

### Y.mat

Y.mat <- read.table(paste0(cca_dir, "data/Y.txt"))
df_y <- Y.mat %>% 
    rownames_to_column("sample") %>% 
    as_tibble()
covariates <- df_y %>% dplyr::select(-sample) %>% colnames

### X.mat

X.mat <- read.table(paste0(cca_dir, "data/X.txt"))
df_x <- X.mat %>% 
    t() %>% 
    as.data.frame %>% 
    rownames_to_column("transcript_id") %>% 
    as_tibble()
transcripts <- df_x %>% pull(transcript_id) 
length(transcripts) # n = 54302 (_WITHGRAY); n = 24626 (_withGray)

n_gray <- df_modules_filt %>% filter(module == 0) %>% nrow
nrow(df_modules_filt) - n_gray


# model weights & standard deviations -------------------------------------

### Y
df_y_weights <- model_1$wY %>% 
    as.data.frame %>% 
    as_tibble() %>% 
    .[best_split,] %>% 
    pivot_longer(1:ncol(.), names_to = "covariate", values_to = "weight") %>% 
    mutate(covariate = covariates) %>% 
    
    ## IMPORTANT: negate weights so SCZ is positive (easier to understand)
    mutate(weight = -weight)

df_y_sd <- boot_1$wY[best_split, , ] %>% 
    as.data.frame %>% 
    as_tibble() %>% 
    rename_all(~ covariates) %>% 
    map_dfr( ~ sd(.x)) %>% 
    pivot_longer(1:ncol(.), names_to = "covariate", values_to = "sd")

df_y_weight_sd <- df_y_weights %>% 
    left_join(df_y_sd) %>% 
    mutate(z_score = weight/sd) %>% 
    mutate(covariate = case_when(
        covariate == "control" ~ "Control",
        covariate == "bipolar" ~ "BD",
        covariate == "mdd" ~ "MDD",
        covariate == "schizo" ~ "SCZ",
        str_detect(covariate, "dim") ~ str_replace(covariate, "dim", "Dim "),
        TRUE ~ covariate
    )
    ) %>% 
    arrange(-abs(z_score))

### X
df_x_weights <- model_1$wX %>% 
    as.data.frame %>% 
    as_tibble() %>% 
    .[best_split,] %>% 
    pivot_longer(1:ncol(.), names_to = "ensembl_transcript_id", values_to = "weight") %>% 
    mutate(ensembl_transcript_id = transcripts) %>% 
    
    ## IMPORTANT: negate weights so SCZ is positive (easier to understand)
    mutate(weight = -weight)

df_x_sd <- boot_1$wX[best_split, , ] %>% 
    as.data.frame %>% 
    as_tibble() %>% 
    rename_all(~ transcripts) %>% 
    map_dfr( ~ sd(.x)) %>% 
    pivot_longer(1:ncol(.), names_to = "ensembl_transcript_id", values_to = "sd")

df_x_weight_sd <- df_x_weights %>% 
    left_join(df_x_sd) %>% 
    mutate(z_score = weight/sd) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(transcript_symbol = ifelse(is.na(transcript_symbol), ensembl_transcript_id, transcript_symbol)) %>% 
    left_join(df_modules_filt) %>% 
    dplyr::select(module, color, ensembl_transcript_id, transcript_symbol, weight, sd, z_score) %>% 
    arrange(module, -abs(z_score))


# Calculate latent variables and structural correlations ------------------

### CALCULATE LATENT VARIABLE
df_lvs <- tibble(lvx = (as.matrix(X.mat) %*% (df_x_weights$weight))[,1],
                 lvy = (as.matrix(Y.mat) %*% (df_y_weights$weight))[,1]) %>% 
    mutate(sample = rownames(Y.mat), .before = 1) %>% 
    mutate(dx = case_when(
        str_detect(sample, "control") ~ "Control",
        str_detect(sample, "bipolar") ~ "BD",
        str_detect(sample, "mdd") ~ "MDD",
        str_detect(sample, "schizo") ~ "SCZ"
    ) %>% factor(levels = names(dx_colors))
    )

### CALCULATE STRUCTURE CORRELATIONS TO IDENTIFY COVARIATE-GENE EXPRESSION ASSOCIATION

# Y
df_y_struct_cor <- map_dfr(.x = Y.mat, 
                           .f = ~ tibble(
                               pearsons_r = cor.test(.x, df_lvs$lvy)$estimate,
                               p_value = cor.test(.x, df_lvs$lvy)$p.value
                           )
) %>% 
    mutate(covariate = names(Y.mat), .before = 1) %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr")#,
           #significant = ifelse(p_adj < 0.05, "*", "-")
    ) %>% 
    mutate(covariate = case_when(
        covariate == "control" ~ "Control",
        covariate == "bipolar" ~ "BD",
        covariate == "mdd" ~ "MDD",
        covariate == "schizo" ~ "SCZ",
        str_detect(covariate, "dim") ~ str_replace(covariate, "dim", "Dim "),
        TRUE ~ covariate
    ))  

df_y_res <- df_y_weight_sd %>% 
    left_join(df_y_struct_cor) %>% 
    arrange(p_adj, -abs(z_score)) %>% 
    mutate(significant = ifelse(p_adj < 0.05 & abs(z_score) >= 2, 1, 0))

# X
df_x_struct_cor <- map_dfr(.x = X.mat, 
                           .f = ~ tibble(
                               pearsons_r = cor.test(.x, df_lvs$lvx)$estimate,
                               p_value = cor.test(.x, df_lvs$lvx)$p.value
                           )
) %>% 
    mutate(ensembl_transcript_id = names(X.mat), .before = 1) %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    left_join(df_modules_filt) %>% 
    mutate(transcript_symbol = ifelse(is.na(transcript_symbol), ensembl_transcript_id, transcript_symbol))

df_x_res <- df_x_weight_sd %>% 
    left_join(df_x_struct_cor) %>% 
    dplyr::select(-color) %>% 
    arrange(p_adj, -abs(z_score)) %>% 
    dplyr::select(ensembl_gene_id, gene_symbol, everything()) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    mutate(significant = ifelse(p_adj < 0.05 & abs(z_score) >= 2, 1, 0))



# Export tables for publication -------------------------------------------

# to project dir
# write_xlsx(list("A. Model results across VARx" = as.data.frame(df_results),
#                 "B. Covariate (y) weights" = as.data.frame(df_y_res),
#                 "C. Gene (x) weights" = as.data.frame(df_x_res)
# ),
# paste0(base_dir, "outputs/tables/for_manuscript/TableS7_TRANSCRIPTS_GRCCA_res_withGray.xlsx")
# )

