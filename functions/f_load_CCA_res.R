
########################################################################################

## Function that loads results from the RCCA toolkit and formats result into a table ##

## INPUTS
## - cca_directory = directory that contains CCA data & framework folders
## - type = specify if genes or transcripts (character)
## - analysis_type = CCA algorithm used, either rcca or grcca
## - mu = mu (x group regularization parameter) value used in analysis (numeric or character ok; set to NULL if RCCA was run (though any value in there is fine, it won't be used))
## - VARx = was variance explained in the X matrix optimized ("0.1_1") or set at 0.99?
## - include_Cmat = logical indicating whether Cmat was included in the analysis (were confounds regressed?)

########################################################################################

f_load_cca_res <- function(cca_directory = cca_dir,
                           level = c("gene", "transcript"),
                           analysis_type = c("rcca", "grcca"),
                           mu,
                           lambda = "Nfeat",
                           VARx = c("0.1_1", "0.99"),
                           include_Cmat = c(TRUE, FALSE),
                           rename_covariates = c(TRUE, FALSE),
                           all_effects = c(TRUE, FALSE)
                           ) {
    
### Set paths & data split to load results ###
    
    print("Setting paths")
    
    ## Identify correct framework analysis directory
    if (include_Cmat) {cmat <- "withCmat"} else {cmat <- "noCmat"}
    if (analysis_type == "grcca") {group <- paste0("_groupL2x", as.character(mu))} else {group <- ""}
    if (all_effects == TRUE) {analysis_dir <- paste0(analysis_type, "_permutation_VARx", VARx, "_L2x", lambda, group, "_", cmat, "_allEffects")
    } else {analysis_dir <- paste0(analysis_type, "_permutation_VARx", VARx, "_L2x", lambda, group, "_", cmat)}
    
    ## Combine paths to point to results (and check to make sure directory exists)
    results_dir <- paste0(cca_directory, "framework/", analysis_dir)
    if (!dir.exists(results_dir)) {stop("Analysis directory does not exist\n", results_dir)}
    
    ## Identify final model path (weights)
    model_path <- paste0(results_dir, "/res/level1/model_1.mat")
    
    # Identify bath to bootstrapped weights (to calculate sd)
    boot_path <- paste0(results_dir, "/boot/level1/allboot_1.mat")
    
    ## Load results table, include X-Y LV correlation across X variance explain values
    df_results <- read_table(paste0(results_dir, "/res/level1/results_table.txt"))
    
    ## Select the split (defined by variance explained) with the lowest p-value & highest X-Y latent variable correlation
    best_split <- df_results %>% arrange(pval, -correl) %>% pull(set) %>% .[[1]]
    
    print(paste0("Model optimized at ", best_split*10, "% variance"))
  
### Load X.mat & Y.mat ###
    
    print("Loading Y and X matrices")
    
    ## Y.mat
    Y.mat <- read.table(paste0(cca_directory, "data/Y.txt"))
    df_y <- Y.mat %>% 
        rownames_to_column("sample") %>% 
        as_tibble()
    covariates <- df_y %>% dplyr::select(-sample) %>% colnames
    
    # rename covariates for plotting
    if (rename_covariates == TRUE) {
        covariate_names <- case_when(
            covariates == "control" ~ "Control",
            covariates == "bipolar" ~ "BD",
            covariates == "mdd" ~ "MDD",
            covariates == "schizo" ~ "SCZ",
            str_detect(covariates, "dim") ~ str_replace(covariates, "dim", "Dim ") 
        )
        names(Y.mat) <- covariate_names
    } else if (rename_covariates == FALSE) {
        covariate_names <- names(Y.mat) %>% str_replace("dim", "Dim ")
        names(Y.mat) <- covariate_names
    }
    
    ## X.mat
    if (level == "gene") {col_name <- "ensembl_gene_id"} else if (level == "transcript") {col_name <- "ensembl_transcript_id"}
    X.mat <- read.table(paste0(cca_directory, "data/X.txt"))
    df_x <- X.mat %>% 
        t() %>% 
        as.data.frame %>% 
        rownames_to_column(col_name) %>% 
        as_tibble()
    gene_order <- df_x[[paste0(col_name)]]

### Read & format model weights & standard deviations ###
    
    print("Loading model weights & calculating standard deviations")
    
    ## Y
    
    # weights
    y_weights <- h5read(model_path, name = "wY")[best_split, ]
    names(y_weights) <- covariate_names
    df_y_weights <- enframe(y_weights, name = "covariate", value = "weight")
    
    # bootstrapped weights
    y_boots <- h5read(boot_path, name = "wY")[best_split, , ]
    colnames(y_boots) <- covariate_names
    df_y_sd <- apply(y_boots, 2, sd) %>% 
        enframe(name = "covariate", value = "sd")
    
    # calculate z-score
    df_y_weight_sd <- df_y_weights %>% 
        left_join(df_y_sd) %>% 
        mutate(z_score = weight/sd) %>% 
        arrange(-abs(z_score))
    
    ## X
    
    # weights
    x_weights <- h5read(model_path, name = "wX")[best_split, ]
    names(x_weights) <- gene_order
    df_x_weights <- enframe(x_weights, name = col_name, value = "weight")
    
    # bootstrapped weights
    x_boots <- h5read(boot_path, name = "wX")[best_split, , ]
    colnames(x_boots) <- gene_order
    df_x_sd <- apply(x_boots, 2, sd) %>% 
        enframe(name = col_name, value = "sd")
    
    # calculate z-score
    df_x_weight_sd <- df_x_weights %>% 
        left_join(df_x_sd) %>% 
        mutate(z_score = weight/sd) %>% 
        arrange(-abs(z_score))
  
### Calculate latent variables ###
    
    print("Calculating latent variables")
    
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
  
### Calculate structure correlations to identify covariate-gene expression association ###
    
    print("Calculating structure correlations")
    
    ## Y
    df_y_struct_cor <- map_dfr(.x = Y.mat, 
                               .f = ~ tibble(
                                   pearsons_r = cor.test(.x, df_lvs$lvy)$estimate,
                                   p_value = cor.test(.x, df_lvs$lvy)$p.value
                               )
    ) %>% 
        mutate(covariate = names(Y.mat), .before = 1) %>% 
        mutate(p_adj = p.adjust(p_value, method = "fdr")) 
    
    ## X
    df_x_struct_cor <- map_dfr(.x = X.mat, 
                               .f = ~ tibble(
                                   pearsons_r = cor.test(.x, df_lvs$lvx)$estimate,
                                   p_value = cor.test(.x, df_lvs$lvx)$p.value
                               )
    ) %>% 
        mutate(!!col_name := names(X.mat), .before = 1) %>% 
        mutate(p_adj = p.adjust(p_value, method = "fdr"))
  
### Combine z-score and structure correlations into X & Y results tables ###
    
    print("Generating results tables")
    
    ## Y
    df_y_res <- df_y_weight_sd %>% 
        left_join(df_y_struct_cor) %>% 
        arrange(p_adj, -abs(z_score)) %>% 
        mutate(significant = ifelse(p_adj < 0.05 & abs(z_score) >= 1.96, 1, 0))
    
    ## X
    df_x_res <- df_x_weight_sd %>% 
        left_join(df_x_struct_cor) %>% 
        arrange(p_adj, -abs(z_score)) %>% 
        mutate(significant = ifelse(p_adj < 0.05 & abs(z_score) >= 1.96, 1, 0))
    
 
### Return model information and results tibbles as list ###
    
    return(
        list(
            "mu" = mu,
            "confounds_regressed" = include_Cmat,
            "latent_variables" = df_lvs,
            "model_results" = df_results,
            "y_res" = df_y_res,
            "x_res" = df_x_res
        )
    )
    
}
