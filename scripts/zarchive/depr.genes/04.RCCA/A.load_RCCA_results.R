
########################################################################################

# Load gene-level RCCA results to analyze

########################################################################################


# Setup -------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/load_WGCNA_res.R"))


# Set paths & data split to load results -----------------------------------------------------


## Set paths to read in RCCA res (generated in cluster)
n_mca_dim <- 8
type <- "GENES/"
project_dir <- paste0(prefix,
                      "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "_", 
                      n_mca_dim, "MCA_regressBrainWeight_WITHGRAY/") # _WITHGRAY, _withGray
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)
analysis_dir <- "rcca_permutation_VARx0.1_1_L2xNfeat_withCmat_allEffects"
model_path <- paste0(cca_dir, "framework/", analysis_dir, "/res/level1/model_1.mat")
boot_path <- paste0(cca_dir, "framework/", analysis_dir, "/boot/level1/allboot_1.mat")

## Read results table output
df_results <- read_table(paste0(cca_dir, "framework/", analysis_dir, "/res/level1/results_table.txt"))

## Select the split (defined by variance explained) with the lowest p-value &
## highest X-Y latent variable correlation
best_split <- df_results %>% arrange(pval, -correl) %>% pull(set) %>% .[[1]]


# Load X.mat & Y.mat ------------------------------------------------------

## Y.mat
Y.mat <- read.table(paste0(cca_dir, "data/Y.txt"))
df_y <- Y.mat %>% 
    rownames_to_column("sample") %>% 
    as_tibble()
covariates <- df_y %>% dplyr::select(-sample) %>% colnames

# rename covariates for plotting
covariates_renamed <- case_when(
    covariates == "control" ~ "Control",
    covariates == "bipolar" ~ "BD",
    covariates == "mdd" ~ "MDD",
    covariates == "schizo" ~ "SCZ",
    str_detect(covariates, "dim") ~ str_replace(covariates, "dim", "Dim ") 
)
names(Y.mat) <- covariates_renamed

## X.mat
X.mat <- read.table(paste0(cca_dir, "data/X.txt"))
df_x <- X.mat %>% 
    t() %>% 
    as.data.frame %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    as_tibble()
gene_order <- df_x$ensembl_gene_id
length(gene_order) # n = 18677 _WITHGRAY
all(gene_order == colnames(X.mat))


# Read & format model weights & standard deviations -------------------------------------

### Y ###

# weights
y_weights <- h5read(model_path, name = "wY")[best_split, ]
names(y_weights) <- covariates_renamed
df_y_weights <- enframe(y_weights, name = "covariate", value = "weight") %>% 
    
    ## IMPORTANT: negate weights so SCZ is positive (easier to understand)
    mutate(weight = -weight)

# bootstrapped weights
y_boots <- h5read(boot_path, name = "wY")[best_split, , ]
colnames(y_boots) <- covariates_renamed
df_y_sd <- apply(y_boots, 2, sd) %>% 
    enframe(name = "covariate", value = "sd")

# calculate z-score
df_y_weight_sd <- df_y_weights %>% 
    left_join(df_y_sd) %>% 
    mutate(z_score = weight/sd) %>% 
    arrange(-abs(z_score))

### X ###

# weights
x_weights <- h5read(model_path, name = "wX")[best_split, ]
names(x_weights) <- gene_order
df_x_weights <- enframe(x_weights, name = "ensembl_gene_id", value = "weight") %>% 
    
    ## IMPORTANT: negate weights so SCZ is positive (easier to understand)
    mutate(weight = -weight)

# bootstrapped weights
x_boots <- h5read(boot_path, name = "wX")[best_split, , ]
colnames(x_boots) <- gene_order
df_x_sd <- apply(x_boots, 2, sd) %>% 
    enframe(name = "ensembl_gene_id", value = "sd")

# calculate z-score
df_x_weight_sd <- df_x_weights %>% 
    left_join(df_x_sd) %>% 
    mutate(z_score = weight/sd) %>% 
    arrange(-abs(z_score))



# Calculate latent variables and structural correlations ------------------

## Calculate latent variables
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

## Calculate structure correlations to identify covariate-gene expression association

# Y
df_y_struct_cor <- map_dfr(.x = Y.mat, 
                           .f = ~ tibble(
                               pearsons_r = cor.test(.x, df_lvs$lvy)$estimate,
                               p_value = cor.test(.x, df_lvs$lvy)$p.value
                           )
) %>% 
    mutate(covariate = names(Y.mat), .before = 1) %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr")) 

# X
df_x_struct_cor <- map_dfr(.x = X.mat, 
                           .f = ~ tibble(
                               pearsons_r = cor.test(.x, df_lvs$lvx)$estimate,
                               p_value = cor.test(.x, df_lvs$lvx)$p.value
                           )
) %>% 
    mutate(ensembl_gene_id = names(X.mat), .before = 1) %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr"))


# Combine z-score and structure correlations into results table -----------

## Y
df_y_res <- df_y_weight_sd %>% 
    left_join(df_y_struct_cor) %>% 
    arrange(p_adj, -abs(z_score)) %>% 
    mutate(significant = ifelse(p_adj < 0.05 & abs(z_score) >= 2, 1, 0))

## X
df_x_res <- df_x_weight_sd %>% 
    left_join(df_x_struct_cor) %>% 
    arrange(p_adj, -abs(z_score)) %>% 
    mutate(significant = ifelse(p_adj < 0.05 & abs(z_score) >= 2, 1, 0)) %>% 
    
    # add gene symbol and module information
    left_join(df_ensembl_to_symbol, by = join_by(ensembl_gene_id)) %>% 
    left_join(df_modules_filt %>% dplyr::select(-c(mod_set, color)), by = join_by(ensembl_gene_id)) %>% 
    dplyr::select(ensembl_gene_id, gene_symbol, module, everything()) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol))
    

# Export tables -----------------------------------------------------------


tables_dir <- paste0(base_dir, "outputs/tables/02Aug_update/")
write_xlsx(list("A. Model results across VARx" = as.data.frame(df_results),
                "B. Covariate (y) weights" = as.data.frame(df_y_res),
                "C. Gene (x) weights" = as.data.frame(df_x_res)
),
paste0(tables_dir, "gene_RCCA_results.xlsx")
)

