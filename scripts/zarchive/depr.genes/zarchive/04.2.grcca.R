
########################################################################################

# Explore GRCCA

########################################################################################



# libraries ---------------------------------------------------------------


# remotes::install_github("ElenaTuzhilina/RCCA")

  library(tidyverse)
  library(patchwork)
  library(janitor)
  library(RColorBrewer)
  library(RCCA)
  library(biomaRt)
  library(tidymodels)  
  library(ggrepel)
  library(raveio)
  library(ggpubr)
  library(ggchicklet)



# set theme for plots -----------------------------------------------------


  theme_set(theme_bw() +
              theme(plot.title = element_text(size = 18),
                    axis.title = element_text(size = 15),
                    axis.text = element_text(size = 12),
                    strip.text = element_text(size = 15),
                    legend.title = element_text(size = 15),
                    legend.text = element_text(size = 12)))
  



# data --------------------------------------------------------------------


soft_power <- 3
base_dir <- "/Users/work/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/projects/WGCNA/acsg_wgcna/"
prefix <- "20230228_185samples_19kgenes_vst_qSVA123567_MP_RNAbatch_Race_resids"

# LOAD OBJECTS  
  load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
  load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

# IDENTIFY MODULE SET OF INTEREST  
  df_modules_filt <- df_modules %>% 
    filter(min_size == 40 & cut_height == 0.97) %>% 
    unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module)

# COVARIATES
  load("data/covariates/185_all_covariates_clean.Rdata")

# GENE SYMBOL TO ENSEMBL GENE ID MAPPING
  load("objects/ensembl_id_to_gene_symbol.Rdata") # df_ensembl_to_symbol

  

# explain GRCCA hyperparameters -------------------------------------------

  grcca_hyperparams_dir <- paste0(base_dir, "GRCCA_figs/GRCCA_hyperparams/")
  mu0_lambda0 <- read_mat(paste0(grcca_hyperparams_dir, "framework/grcca_holdout1-0.20_mu0_lambda0/res/level1/model_1.mat"))
  mu0_lambda1 <- read_mat(paste0(grcca_hyperparams_dir, "framework/grcca_holdout1-0.20_mu0_lambda1/res/level1/model_1.mat"))
  mu1_lambda0 <- read_mat(paste0(grcca_hyperparams_dir, "framework/grcca_holdout1-0.20_mu1_lambda0/res/level1/model_1.mat"))
  mu1_lambda1 <- read_mat(paste0(grcca_hyperparams_dir, "framework/grcca_holdout1-0.20_mu1_lambda1/res/level1/model_1.mat"))
  
  df_modules <- read_csv(paste0(grcca_hyperparams_dir, "data/LabelsX.csv")) %>% 
    mutate(module = factor(Category, levels = 0:max(Category)))
  
  df_x_hyperparam <- tibble(ensembl_gene_id = colnames(X.mat),
                            mu0_lambda0 = mu0_lambda0$wX,
                            mu0_lambda1 = mu0_lambda1$wX,
                            mu1_lambda0 = mu1_lambda0$wX,
                            mu1_lambda1 = mu1_lambda1$wX) %>% 
    bind_cols(df_modules) %>% dplyr::select(-Label, -Category) %>% 
    pivot_longer(starts_with("mu"), names_to = "framework", values_to = "weight") %>% 
    arrange(module) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = unique(.$ensembl_gene_id)))
  
# PLOT  
  df_x_hyperparam %>% 
    ggplot(aes(x = ensembl_gene_id, y = weight, fill = module)) +
    geom_col() +
    scale_fill_manual(values = c(df_modules_filt %>% 
                                   dplyr::count(module, color) %>% 
                                   dplyr::select(-n) %>% 
                                   pull(color))
    ) +
    facet_wrap(vars(framework), scales = "free_y") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_blank()) +
    ggtitle("X coefficient weights across hyperparameter combinations")
  
  
  df_x_hyperparam %>% 
    filter(str_detect(framework, "lambda0")) %>% 
    group_by(ensembl_gene_id) %>% 
    mutate(mean_weight = mean(weight)) %>% 
    
    ggplot(aes(x = mean_weight)) +
    geom_histogram()
  
# MU
  tibble(
    
    x = 1:1000,
    group = c(rep("a", 375), 
              rep("b", 250), 
              rep("c", 375)),
    
    y1 = c(rnorm(375, mean = 5, sd = 2), 
           rnorm(250, mean = 2, sd = 2), 
           rnorm(375, mean = -4, sd = 2)),
    
    y2 = c(rnorm(375, mean = 3, sd = 1.5), 
           rnorm(250, mean = 0, sd = 1.5), 
           rnorm(375, mean = -2, sd = 1.5)),
    
  ) %>% 
    pivot_longer(starts_with("y"), names_to = "facet", values_to = "y") %>% 
    
  ggplot(aes(x = x, y = y)) +
    geom_col(aes(alpha = y, color = group)) +
    scale_color_manual(values = c("#ffda80ff", "#80ffd7ff", "#fa86f8f5")) +
    
    facet_wrap(vars(facet)) +
    labs(x = "", y = "weights") +
    #ggtitle("shrinkage regularization") +
    theme(strip.text = element_blank(),
          axis.text = element_blank(),
          legend.position = "none")
  
# LAMBDA
  tibble(
    
    x = 1:1000,
    group = c(rep("a", 375), 
              rep("b", 250), 
              rep("c", 375)),
    
    y1 = c(rnorm(375, mean = 5, sd = 2), 
           rnorm(250, mean = 2, sd = 2), 
           rnorm(375, mean = -4, sd = 2)),
    
    y2 = c(rnorm(375, mean = 5, sd = 0.5), 
           rnorm(250, mean = 2, sd = 0.5), 
           rnorm(375, mean = -4, sd = 0.5)),
    
  ) %>% 
    pivot_longer(starts_with("y"), names_to = "facet", values_to = "y") %>% 
    
    ggplot(aes(x = x, y = y)) +
    geom_col(aes(alpha = y, color = group)) +
    scale_color_manual(values = c("#ffda80ff", "#80ffd7ff", "#fa86f8f5")) +
    
    facet_wrap(vars(facet)) +
    labs(x = "", y = "weights") +
    #ggtitle("shrinkage regularization") +
    theme(strip.text = element_blank(),
          axis.text = element_blank(),
          legend.position = "none")
  
  
  
# define matrices ---------------------------------------------------------------

  
# Select covariates based on significant covariates from 4.1: 
  # anti-anxiety, opioids, brain_weight, age_death
  
# DEFINE MATRICES  
  df_grcca_x <- df_vsd_regress %>% 
    dplyr::select(ensembl_gene_id, sample, resids) %>% 
    left_join(df_modules_filt) %>% 
    arrange(module) 
  
  X.mat <- df_grcca_x %>% 
    pivot_wider(id_cols = sample, names_from = ensembl_gene_id, values_from = resids) %>% 
    as.data.frame %>% 
    column_to_rownames("sample")
  
  x_group <- df_grcca_x %>% 
    dplyr::select(ensembl_gene_id, module) %>% 
    distinct() %>% 
    pull(module) %>% 
    paste0("M", .)
  
  # filter x for no grey module genes
  #X.mat <- X.mat[x_group != "M0"]
  #x_group <- x_group[x_group != "M0"]
  
  # remove M for Agoston's toolbox
  x_group <- str_remove(x_group, "M") %>% as.numeric
  
  YC.mat <- df_covariates_all_clean %>% 
    dplyr::select(sample, opioids, anti_anxiety, brain_weight, age_death) %>%  
    mutate(control = ifelse(grepl("control", sample), 1, 0),
           bipolar = ifelse(grepl("bipolar", sample), 1, 0),
           mdd = ifelse(grepl("mdd", sample), 1, 0),
           schizo = ifelse(grepl("schizo", sample), 1, 0),
           
           anti_anxiety = ifelse(anti_anxiety == "positive", 1, 0), 
           opioids = ifelse(opioids == "positive", 1, 0)#,
           ) %>% 
    dplyr::select(sample, control, bipolar, mdd, schizo, 
                  opioids, anti_anxiety, brain_weight, age_death) %>%
    column_to_rownames("sample") %>% 
    as.data.frame
  
  Y.mat <- YC.mat %>% dplyr::select(-c(control, brain_weight, age_death)) # remove control column so columns aren't linearly equivalent
  
  C.mat <- YC.mat %>% dplyr::select(-c(control, bipolar, mdd, schizo, opioids, anti_anxiety))
  
# EXPORT FOR AGOSTON'S TOOLBOX
  # cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_19kgenes_covarOpioidsAnx_regressBrainweightAgedeath/")
  # dir.create(cca_dir)
  # cca_data_dir <- paste0(cca_dir, "data/")
  # dir.create(cca_data_dir)
  # 
  # write.table(X.mat, file = paste0(cca_data_dir, "X.txt"))
  # write.table(Y.mat, file = paste0(cca_data_dir, "Y.txt"))
  # write.table(C.mat, file = paste0(cca_data_dir, "C.txt"))
  # write.table(x_group, file = paste0(cca_data_dir, "XGroup.txt"))
  
  

# GRCCA results -----------------------------------------------------------

  # load("objects/df_hsapiens_genome.RDS") # df_hsapiens_genome
  # 
  # df_ensembl_to_symbol <- df_hsapiens_genome %>% 
  #   dplyr::select(gene_id, gene_name) %>% 
  #   distinct() %>% 
  #   dplyr::rename("ensembl_gene_id" = "gene_id",
  #                 "gene_symbol" = "gene_name") %>% 
  #   mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol))
  # 
  # save(df_ensembl_to_symbol, file = "objects/ensembl_id_to_gene_symbol.Rdata")
  
# GRCCA RES  
  cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_19kgenes_covarOpioidsAnx_regressBrainweightAgedeath/")
  model_mat1 <- read_mat(paste0(cca_dir, "framework/grcca_holdout1-0.20_subsamp5-0.20/res/level1/model_1.mat"))
  
  
# LEVEL 1 X-Y COR
  df_lvs1 <- tibble(lvx = as.matrix(X.mat) %*% model_mat1$wX,
                    lvy = as.matrix(Y.mat) %*% model_mat1$wY)
  
  df_lvs1 %>% 
    ggplot(aes(x = lvx, y = lvy)) +
    geom_point() +
    geom_smooth(method = lm, lty = 2) +
    xlab("X latent variable (x matrix • x weights)") +
    ylab("Y latent variable (y matrix • y weights)") +
    #stat_regline_equation() +
    stat_cor(label.y = 0.00, label.x = -3) +
    ggtitle("Gene x-y latent variable correlation")
  
  
# LEVEL 1 Y   
  df_y_res1 <- tibble(covariate = colnames(Y.mat), 
                      weight = model_mat1$wY) %>% 
    mutate(covariate = factor(covariate, levels = c("bipolar", "mdd", "schizo", 
                                                    "opioids", "anti_anxiety")))
  
  df_y_res1 %>% 
    ggplot(aes(x = covariate, y = weight)) +
    geom_point(size = 5, aes(color = weight)) +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    scale_color_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                          limits = c(-0.1, 0.1)) +
    ylim(c(-0.1, 0.1)) +
    xlab("") +
    ggtitle("Y coefficient weights") +
    theme(legend.position = "none")
  
# LEVEL 1 X  
  df_x_res1 <- tibble(ensembl_gene_id = colnames(X.mat),
                      weight = model_mat1$wX) %>% 
    left_join(df_modules_filt) %>% 
    mutate(abs_weight = abs(weight)) %>% 
    arrange(-abs_weight) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = unique(.$ensembl_gene_id)))
  
  df_x_res1 %>% 
    #filter(abs_weight > 0.0005) %>%
    ggplot(aes(x = ensembl_gene_id, y = abs_weight, fill = module)) +
    geom_col() +
    scale_fill_manual(values = c(df_x_res1 %>% 
                                    #head(100) %>% 
                                    dplyr::count(module, color) %>% 
                                    dplyr::select(-n) %>% 
                                    pull(color))
    ) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("absolute value of X coefficient weights")
  
  df_x_res1 %>% 
    slice_max(abs_weight, n = 50) %>% 
    arrange(weight) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = .$gene_symbol)) %>% 
    ggplot(aes(x = weight, y = gene_symbol, fill = weight)) +
    geom_col() +
    #coord_flip() +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                         limits = c(-0.016, 0.015)) +
    ylab("") +
    theme_classic() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12)
    ) +
    ggtitle("Top 20 absolute gene weights")

# MODULE WEIGHTS
  df_x_res1 %>% 
    group_by(module) %>% 
    summarise(mean_weight = mean(abs_weight)) %>% 
    arrange(-mean_weight) %>% 
    print(n = nrow(.))
  
# COMPARE TO SCZ GENOME  
  df_scz_genome <- readxl::read_excel("data/eva/lago_bahn_2022_scz_genome.xlsx") %>% 
    row_to_names(1) %>% 
    clean_names
  lago_bahn_risk_genes <- df_scz_genome$risk_gene_hgnc %>% unique()
  
  gene_universe <- df_x_res1$gene_symbol
  
  lago_bahn_risk_genes_in_universe <- intersect(gene_universe, lago_bahn_risk_genes) 
  
  top_cca_genes <- df_x_res1 %>%
    group_by(module) %>% 
    slice_max(abs_weight, n = 20) %>%
    filter(module %in% c(9, 15, 12)) %>% 
    pull(gene_symbol) # significant
  
  overlap_genes <- intersect(top_cca_genes, lago_bahn_risk_genes_in_universe)
  
  q <- length(overlap_genes) - 1
  m <- length(lago_bahn_risk_genes_in_universe)
  n <- length(gene_universe) - m
  k <- length(top_cca_genes)
  phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  
  df_hyper_res <- df_x_res1 %>% 
    expand_grid(threshold = c(0.01, 0.02, 0.03)) %>% 
    group_by(module, threshold) %>% 
    nest() %>% 
    
    mutate(top_genes = map(
      
      .x = data,
      .y = threshold,
      .f = ~ .x %>% 
        slice_max(abs_weight, n = floor(.y*nrow(.x))) %>% 
        pull(gene_symbol)
      
    ),
    
    overlap_genes = map(
      
      .x = top_genes,
      .f = ~ intersect(.x, lago_bahn_risk_genes_in_universe)
      
    ),
    
    hyper_p_val = map(
      
      .x = top_genes,
      .y = overlap_genes,
      .f = ~ phyper(length(.y) - 1,
                  length(lago_bahn_risk_genes_in_universe),
                  length(gene_universe) - length(lago_bahn_risk_genes_in_universe),
                  length(.x),
                  lower.tail = FALSE, log.p = FALSE)
      
    )
    
    )
  
  df_hyper_res %>% 
    unnest(cols = c(hyper_p_val)) %>% 
    filter(hyper_p_val < 0.05)
  
  # plot overlap genes
  df_x_res1 %>% 
    filter(gene_symbol %in% overlap_genes) %>% 
    arrange(weight) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = .$gene_symbol)) %>% 
    ggplot(aes(x = weight, y = gene_symbol, fill = weight)) +
    geom_col() +
    #coord_flip() +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                         limits = c(-0.012, 0.01)) +
    ylab("") +
    theme_classic() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12)
    ) +
    ggtitle("SCZ druggable targets; top 1% CCA-PCA genes")
  
# SCZ GWAS?
  trubetskoy_genes <- readxl::read_excel("data/eva/trubetskoy_et_al_2022_scz_GWAS.xls") %>% 
    clean_names %>%
    dplyr::select(genes_all) %>% 
    separate(genes_all, into = paste0('gene', seq_len(max(str_count(.$genes_all, ',')+1))), 
             sep = ",", fill = "right") %>% 
    pivot_longer(1:ncol(.), names_to = "gene", values_to = "gene_symbol") %>% 
    filter(!is.na(gene_symbol)) %>% 
    pull(gene_symbol) %>% 
    unique()
    
  trubetskoy_genes_in_universe <- intersect(gene_universe, trubetskoy_genes) 
  
  top_cca_genes <- df_x_res1 %>%
    group_by(module) %>% 
    slice_max(abs_weight, n = 20) %>%
    filter(module %in% c(9, 15, 12, 21, 6, 1)) %>% 
    pull(gene_symbol) # significant
  
  overlap_genes <- intersect(top_cca_genes, trubetskoy_genes_in_universe)
  
  q <- length(overlap_genes) - 1
  m <- length(trubetskoy_genes_in_universe)
  n <- length(gene_universe) - m
  k <- length(top_cca_genes)
  phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  
# TRANS-DIAGNOSTIC  
  df_pgc_gwas <- readxl::read_excel("data/PGC_transdiagnostic_GWAS.xlsx") %>% 
    row_to_names(2) %>% 
    clean_names
  
  pgc_genes <- df_pgc_gwas %>% 
    dplyr::select(candidate_gene) %>% 
    mutate(candidate_gene = gsub("\\s*\\([^\\)]+\\)", "", candidate_gene)) %>% 
    filter(candidate_gene != "-") %>% 
    separate(candidate_gene, into = paste0('gene', seq_len(max(str_count(.$candidate_gene, ';')+1))), 
             sep = ";", fill = "right") %>% 
    pivot_longer(1:ncol(.), names_to = "gene", values_to = "gene_symbol") %>% 
    filter(!is.na(gene_symbol)) %>% 
    pull(gene_symbol) %>% 
    unique()
  
  pgc_genes_in_universe <- intersect(gene_universe, pgc_genes) 
  
  # overall?
  top_cca_genes <- df_x_res1 %>%
    slice_max(abs_weight, n = 1000) %>% 
    pull(gene_symbol)
  
  overlap_genes <- intersect(top_cca_genes, pgc_genes_in_universe)
  
  q <- length(overlap_genes) - 1
  m <- length(pgc_genes_in_universe)
  n <- length(gene_universe) - m
  k <- length(top_cca_genes)
  phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  
  # plot positions in weights
  df_x_res1 %>% 
    mutate(gwas = ifelse(gene_symbol %in% pgc_genes, 1, 0) %>% factor) %>% 
    
    ggplot(aes(x = abs_weight, y = ensembl_gene_id)) +
    geom_col(aes(fill = gwas)) +
    scale_fill_manual(values = c("grey", "black")) +
    theme(axis.text = element_blank())
  
  # different thresholds?
  df_hyper_res <- df_x_res1 %>% 
    expand_grid(threshold = c(0.05, 0.10, 0.15)) %>% 
    group_by(module, threshold) %>% 
    nest() %>% 
    
    mutate(top_genes = map(
      
      .x = data,
      .y = threshold,
      .f = ~ .x %>% 
        slice_max(abs_weight, n = floor(.y*nrow(.x))) %>% 
        pull(gene_symbol)
      
    ),
    
    overlap_genes = map(
      
      .x = top_genes,
      .f = ~ intersect(.x, pgc_genes_in_universe)
      
    ),
    
    hyper_p_val = map(
      
      .x = top_genes,
      .y = overlap_genes,
      .f = ~ phyper(length(.y) - 1,
                    length(pgc_genes_in_universe),
                    length(gene_universe) - length(pgc_genes_in_universe),
                    length(.x),
                    lower.tail = FALSE, log.p = FALSE)
      
    )
    
    )
  
  df_hyper_res %>% 
    unnest(cols = c(hyper_p_val)) %>% 
    filter(hyper_p_val < 0.05)
  
  df_hyper_res %>% print(n = nrow(.))
  
  
# CHECK DIFFERENTIAL EXPRESSION
  load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res
  
  top_cca_genes <- df_x_res1 %>%
    head(0.01*nrow(.)) %>%
    pull(ensembl_gene_id) # significant
  
  df_akula_res %>% 
    filter(id %in% top_cca_genes & padj < 0.05)
  
  
  
# PCA-CCA RESULTS -----------------------------------------------


    
  #cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_19kgenes_regressOpioidsSuicideBrainweightAgedeath/")
  #cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_19kgenes_CovarOpioidsSuicide_regressBrainweightAgedeath/")
  #cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_19kgenes_regressOpioidsAnxBrainweightAgedeath/")
  cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_19kgenes_covarOpioidsAnx_regressBrainweightAgedeath/")
  model_mat1 <- read_mat(paste0(cca_dir, "framework/cca_holdout1-0.20_subsamp5-0.20/res/level1/model_1.mat"))
  

# LEVEL 1 X-Y COR
  df_lvs1 <- tibble(lvx = as.matrix(X.mat) %*% model_mat1$wX,
                   lvy = as.matrix(Y.mat) %*% model_mat1$wY)
  
  df_lvs1 %>% 
    ggplot(aes(x = lvx, y = lvy)) +
    geom_point() +
    geom_smooth(method = lm, lty = 2) +
    xlab("X latent variable (x matrix • x weights)") +
    ylab("Y latent variable (y matrix • y weights)") +
    #stat_regline_equation() +
    stat_cor(label.y = 0.13, label.x = -0.2) +
    ggtitle("Gene x-y latent variable correlation")
  
  
# LEVEL 1 Y   
  df_y_res1 <- tibble(covariate = colnames(Y.mat), 
                     weight = model_mat1$wY) %>% 
    mutate(covariate = factor(covariate, levels = c("bipolar", "mdd", "schizo", 
                                                    "opioids", "anti_anxiety")))
  
  df_y_res1 %>% 
    ggplot(aes(x = covariate, y = weight)) +
    geom_point(size = 5, aes(color = weight)) +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    scale_color_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                          limits = c(-0.1, 0.1)) +
    ylim(c(-0.1, 0.1)) +
    xlab("") +
    ggtitle("Y coefficient weights") +
    theme(legend.position = "none")
  
# LEVEL 1 X  
  df_x_res1 <- tibble(ensembl_gene_id = colnames(X.mat),
                     weight = model_mat1$wX) %>% 
    left_join(df_modules_filt) %>% 
    mutate(abs_weight = abs(weight)) %>% 
    arrange(-abs_weight) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = unique(.$ensembl_gene_id)))
    
  df_x_res1 %>% 
    #filter(abs_weight > 0.0005) %>%
    ggplot(aes(x = ensembl_gene_id, y = abs_weight, color = module)) +
    geom_point() +
    scale_color_manual(values = c(df_x_res1 %>% 
                                   #head(100) %>% 
                                   dplyr::count(module, color) %>% 
                                   dplyr::select(-n) %>% 
                                   pull(color))
    ) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("absolute value of X coefficient weights")
    
  df_x_res1 %>% 
    slice_max(abs_weight, n = 20) %>% 
    arrange(weight) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = .$gene_symbol)) %>% 
    ggplot(aes(x = weight, y = gene_symbol, fill = weight)) +
    geom_col() +
    #coord_flip() +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                          limits = c(-0.004, 0.004)) +
    ylab("") +
    theme_classic() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12)
          ) +
    ggtitle("Top 20 absolute gene weights")
  
  df_x_res1 %>% 
    filter(abs_weight > 0.0025) %>% 
    count(module) %>% 
    left_join(df_modules_filt %>% 
                count(module) %>% 
                dplyr::rename("module_n" = "n")
              ) %>% 
    mutate(perc_mod = n/module_n * 100) %>% 
    arrange(-perc_mod)
  

  

# overlap with published results? -----------------------------------------

    
    
# DRUGGABLE TARGETS IN SCZ
    df_scz_genome <- readxl::read_excel("data/eva/lago_bahn_2022_scz_genome.xlsx") %>% 
      row_to_names(1) %>% 
      clean_names
    lago_bahn_risk_genes <- df_scz_genome$risk_gene_hgnc
    
    gene_universe <- df_x_res1$gene_symbol
    
    lago_bahn_risk_genes_in_universe <- intersect(gene_universe, lago_bahn_risk_genes)
    
    top_cca_genes <- df_x_res1 %>%
      #filter(abs_weight > 0.005) %>% 
      head(0.01*nrow(.)) %>%
      pull(gene_symbol) # significant

    overlap_genes <- intersect(top_cca_genes, lago_bahn_risk_genes_in_universe)
    
    q <- length(overlap_genes) - 1
    m <- length(lago_bahn_risk_genes_in_universe)
    n <- length(gene_universe) - m
    k <- length(top_cca_genes)
    
    phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
    
    df_scz_genome %>% filter(risk_gene_hgnc %in% overlap_genes) %>% View()

    

  # plot overlap genes
    df_x_res1 %>% 
     filter(gene_symbol %in% overlap_genes) %>% 
      arrange(weight) %>% 
      mutate(gene_symbol = factor(gene_symbol, levels = .$gene_symbol)) %>% 
      ggplot(aes(x = weight, y = gene_symbol, fill = weight)) +
      geom_col() +
      #coord_flip() +
      scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                           limits = c(-0.003, 0.003)) +
      ylab("") +
      theme_classic() +
      theme(plot.title = element_text(size = 18),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 12),
            strip.text = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12)
      ) +
      ggtitle("SCZ druggable targets; top 1% CCA-PCA genes")
    
# chromosome region
    
    # get gene info from ensembl database
    hsapiens_ensembl <- useEnsembl(biomart = "genes",
                                   dataset = "hsapiens_gene_ensembl")
    
    listAttributes(hsapiens_ensembl) %>% 
      tibble %>% 
      filter(grepl("chrom", name))
    
    df_chromosomes <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                      filters = "hgnc_symbol",
                      values = overlap_genes,
                      mart = hsapiens_ensembl) %>% 
      tibble() %>% 
      filter(chromosome_name != "CHR_HSCHR16_3_CTG1") %>% 
      mutate(chromosome_name = as.numeric(chromosome_name)) %>% 
      arrange(chromosome_name, start_position) %>% 
      left_join(read_csv("data/chromosome_lengths.csv")) %>% 
      mutate(chromosome_name = factor(paste0("chr", chromosome_name), 
                                      levels = paste0("chr", unique(.$chromosome_name)))) 
    
   
    set.seed(123)
    df_chromosomes %>% 

      mutate(xmin = match(chromosome_name, levels(.$chromosome_name)) - 0.3,
             xmax = match(chromosome_name, levels(.$chromosome_name)) + 0.3) %>%
      #mutate(length = ifelse(length > 151000000, 151000000, length)) %>% 

      ggplot(aes(x = chromosome_name)) +
      geom_chicklet(data = df_chromosomes %>% dplyr::select(chromosome_name, length) %>% distinct(),
               mapping = aes(y = length), 
               radius = grid::unit(2, "mm"), width = 0.5,
               color = "black", fill = "white") +
      geom_rect(aes(xmin = xmin, xmax = xmax,
                    ymin = start_position, ymax = end_position,
                    fill = hgnc_symbol)) +
      geom_label_repel(aes(y = start_position, label = hgnc_symbol, fill = hgnc_symbol),
                       alpha = 0.8, min.segment.length = 0.1) +
      scale_fill_manual(values = c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))) +
      #facet_wrap(vars(chromosome_name), scales = "free") +
      coord_flip() +
      labs(x = "", y = "") +
      theme_void() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 12),
            legend.position = "none") +
      ggtitle("Location in genome")
    
    genes_16p11.2 <- c("SPN", "QPRT", "C16ORF54", "ZG16", "KIF22", "MAZ", "PRRT2", "C16orf53",
      "MVP", "CDIPT", "SEZ62L", "ASPHD1", "KCTD13","TMEM219", "TAOK2", "HIRIP3",
      "INO80E", "DOC2A", "C16orf92", "FAM57B", "ALDOA", "PPP4C", "TBX6", "YPEL3", 
      "GDPD3", "MAPK3")

    df_x_res1 %>% 
      filter(gene_symbol %in% genes_16p11.2)
    
    genes_22p11.2 <- c("DGCR2", "DGCR6", "DGCR14", "PRODH", "TSSK2", "GSC2",
                       "SLC25A1", "CLTCL1", "HIRA", "MRPL40", "C22orf39", "UFD1L",
                       "CDC45", "CLDN5", "SEPT5", "GP1BB", "TBX1", "GNB1L", "C22orf29",
                       "TXNRD2", "ARVCF", "TANGO2", "RANBP1", "ZDHHC8", "CCDC188",
                       "LINC00896", "COMT", "DGCR8", "TRMT2A", "RTN4R", "DGCR6L", 
                       "ZNF74", "SCARF2", "KLHL22", "MED15", "CRKL", "SNAP29",
                       "PI4KA", "AIFM3", "SERPIND1", "LZTR1", "THAP7", "P2RX6",
                       "SLC7A4", "HIC2")
    
    df_x_res1 %>% 
      filter(gene_symbol %in% genes_22p11.2)
    
    
    

    
    

# RCCA results ------------------------------------------------------------


################## LEVEL 1    
    
# LEVEL 1 RES    
  cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_19kgenes_covarOpioidsAnx_regressBrainweightAgedeath/")
  model_mat_rcca1 <- read_mat(paste0(cca_dir, "framework/rcca_holdout1-0.20_subsamp5-0.20/res/level1/model_1.mat"))
    
# LEVEL 1 X-Y COR
  df_lvs_rcca1 <- tibble(lvx = as.matrix(X.mat) %*% model_mat_rcca1$wX,
                    lvy = as.matrix(Y.mat) %*% model_mat_rcca1$wY)
  
  df_lvs_rcca1 %>% 
    ggplot(aes(x = lvx, y = lvy)) +
    geom_point(shape = 1) +
    geom_smooth(method = lm, lty = 2, color = "#5961c4", alpha = 0.8) +
    xlab("X latent variable (x matrix • x weights)") +
    ylab("Y latent variable (y matrix • y weights)") +
    #stat_regline_equation() +
    stat_cor(label.y = 0.11, label.x = -3) +
    ggtitle("Gene x-y latent variable correlation")
  
# LEVEL 1 Y   
  df_y_res_rcca1 <- tibble(covariate = colnames(Y.mat), 
                           weight = model_mat_rcca1$wY) %>% 
    mutate(covariate = factor(covariate, levels = c("bipolar", "mdd", "schizo", 
                                                    "opioids", "anti_anxiety")))
  
  df_y_res_rcca1 %>% 
    ggplot(aes(x = covariate, y = weight)) +
    geom_point(size = 5, aes(color = weight)) +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    scale_color_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                          limits = c(-0.1, 0.1)) +
    ylim(c(-0.1, 0.1)) +
    xlab("") +
    ggtitle("Y coefficient weights") +
    theme(legend.position = "none")
  
# LEVEL 1 X  
  df_x_res_rcca1 <- tibble(ensembl_gene_id = colnames(X.mat),
                           weight = model_mat_rcca1$wX) %>% 
    left_join(df_modules_filt) %>% 
    mutate(abs_weight = abs(weight)) %>% 
    arrange(-abs_weight) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) #%>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = unique(.$ensembl_gene_id)))
  
  df_col <- df_modules_filt %>% 
    dplyr::select(module, color) %>% 
    distinct()
  colors <- df_col %>% pull(color)
  names(colors) <- df_col %>% pull(module)
  
  df_x_res_rcca1 %>% 
    arrange(module) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = .$ensembl_gene_id)) %>% 
    
    ggplot(aes(x = ensembl_gene_id, y = weight, fill = module)) +
    geom_col() +
    scale_fill_manual(values = colors) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("absolute value of X coefficient weights")
  
  # plot top n
  df_col <- df_modules_filt %>% 
    dplyr::select(module, color) %>% 
    distinct()
  colors <- df_col %>% pull(color)
  names(colors) <- df_col %>% pull(module)
  
  df_x_res1 %>% 
    slice_max(abs_weight, n = 30) %>% 
    arrange(weight) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = .$gene_symbol)) %>% 
    
    ggplot(aes(x = weight, y = gene_symbol, 
               fill = weight)) +
    geom_col() +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                         limits = c(-0.0151, 0.0151)) +
    ylab("") +
    theme_classic() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12)
    ) +
    ggtitle("Top 30 absolute gene weights")
  
  df_x_res_rcca1 %>% filter(module != 0) %>% 
    head(0.1*nrow(.)) %>% 
    count(module) %>% 
    left_join(df_modules_filt %>% 
                count(module) %>% 
                dplyr::rename("module_n" = "n")
    ) %>% 
    mutate(perc_mod = n/module_n * 100) %>% 
    arrange(-perc_mod)
  
  
# DRUGGABLE TARGETS IN SCZ
  
  # data
  load("objects/df_hsapiens_genome.RDS") # df_hsapiens_genome
  df_scz_genome <- readxl::read_excel("data/eva/lago_bahn_2022_scz_genome.xlsx") %>% 
    row_to_names(1) %>% 
    clean_names
  
  # hypergeometric test
  lago_bahn_risk_genes <- df_scz_genome$risk_gene_hgnc %>% unique()
  
  gene_universe <- df_hsapiens_genome %>% 
    filter(gene_id %in% colnames(X.mat)) %>% 
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>% 
    pull(gene_name) %>% unique
  
  lago_bahn_risk_genes_in_universe <- intersect(gene_universe, lago_bahn_risk_genes) %>% unique
  
  top_genes_rcca1 <- df_x_res_rcca1 %>%
    #filter(abs_weight > 0.005) %>% 
    slice_max(abs_weight, n = floor(0.01*nrow(.))) %>%
    pull(gene_symbol) # significant
  
  overlap_genes1 <- intersect(top_genes_rcca1, lago_bahn_risk_genes_in_universe)
  
  q <- length(overlap_genes1) - 1
  m <- length(lago_bahn_risk_genes_in_universe)
  n <- length(gene_universe) - m
  k <- length(top_genes_rcca1)
  
  phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  
  df_scz_genome %>% filter(risk_gene_hgnc %in% overlap_genes1) %>% View()
  
  # plot overlap genes
  df_x_res_rcca1 %>% 
    filter(gene_symbol %in% overlap_genes1 & abs_weight > 0.005) %>% 
    arrange(weight) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = unique(.$gene_symbol))) %>% 
    ggplot(aes(x = weight, y = gene_symbol, fill = weight)) +
    geom_col() +
    #coord_flip() +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                         limits = c(-0.02, 0.02)) +
    ylab("") +
    theme_classic() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 15),
          legend.position = "none"
    ) +
    ggtitle("SCZ druggable targets; top 1% genes")
  
  # chromosome region
  
  # get gene info from ensembl database
  hsapiens_ensembl <- useEnsembl(biomart = "genes",
                                 dataset = "hsapiens_gene_ensembl")
  
  listAttributes(hsapiens_ensembl) %>% 
    tibble %>% 
    filter(grepl("chrom", name))
  
  df_chromosomes <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                          filters = "hgnc_symbol",
                          values = overlap_genes,
                          mart = hsapiens_ensembl) %>% 
    tibble() %>% 
    mutate(chromosome_name = as.numeric(chromosome_name)) %>% 
    arrange(chromosome_name, start_position) %>% 
    left_join(read_csv("data/chromosome_lengths.csv")) %>% 
    mutate(chromosome_name = factor(paste0("chr", chromosome_name), 
                                    levels = paste0("chr", unique(.$chromosome_name)))) 
  
  
  set.seed(123)
  df_chromosomes %>% 
    left_join(df_modules_filt %>% 
                dplyr::rename("gene_id" = "ensembl_gene_id") %>% 
                left_join(df_hsapiens_genome %>% 
                            dplyr::select(gene_id, gene_name) %>% 
                            distinct() %>% 
                            filter(gene_name %in% overlap_genes)) %>% 
                dplyr::rename("hgnc_symbol" = "gene_name")) %>% 
    
    mutate(xmin = match(chromosome_name, levels(.$chromosome_name)) - 0.3,
           xmax = match(chromosome_name, levels(.$chromosome_name)) + 0.3) %>%
    #mutate(length = ifelse(length > 151000000, 151000000, length)) %>% 
    
    ggplot(aes(x = chromosome_name)) +
    geom_chicklet(data = df_chromosomes %>% 
                    dplyr::select(chromosome_name, length) %>% 
                    distinct(),
                  mapping = aes(y = length), 
                  radius = grid::unit(2, "mm"), width = 0.5,
                  color = "black", fill = "white") +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = start_position, ymax = end_position,
                  fill = hgnc_symbol)) +
    geom_label_repel(aes(y = start_position, label = hgnc_symbol, fill = module),
                     alpha = 0.8, min.segment.length = 0.1) +
    scale_fill_manual(values = colors) +
    #facet_wrap(vars(chromosome_name), scales = "free") +
    coord_flip() +
    labs(x = "", y = "") +
    theme_void() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 12)) +
    ggtitle("Location in genome")
  
  
################## LEVEL 2    
  
# LEVEL 2 RES    
  model_mat_rcca2 <- read_mat(paste0(cca_dir, "framework/rcca_holdout1-0.20_subsamp5-0.20/res/level2/model_1.mat"))
  
# LEVEL 2 X-Y COR
  df_lvs_rcca2 <- tibble(lvx = as.matrix(X.mat) %*% model_mat_rcca2$wX,
                         lvy = as.matrix(Y.mat) %*% model_mat_rcca2$wY)
  
  df_lvs_rcca2 %>% 
    ggplot(aes(x = lvx, y = lvy)) +
    geom_point(shape = 1) +
    geom_smooth(method = lm, lty = 2, color = "#5961c4", alpha = 0.8) +
    xlab("X latent variable (x matrix • x weights)") +
    ylab("Y latent variable (y matrix • y weights)") +
    #stat_regline_equation() +
    stat_cor(label.y = 0.06, label.x = -12) +
    ggtitle("Gene x-y latent variable correlation")
  
  
# LEVEL 2 Y   
  df_y_res_rcca2 <- tibble(covariate = colnames(Y.mat), 
                           weight = model_mat_rcca2$wY) %>% 
    mutate(covariate = factor(covariate, levels = c("bipolar", "mdd", "schizo", 
                                                    "opioids", "anti_anxiety")))
  
  df_y_res_rcca2 %>% 
    ggplot(aes(x = covariate, y = weight)) +
    geom_point(size = 5, aes(color = weight)) +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    scale_color_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                          limits = c(-0.1, 0.1)) +
    ylim(c(-0.1, 0.1)) +
    xlab("") +
    ggtitle("Y coefficient weights") +
    theme(legend.position = "none")
  
# LEVEL 2 X  
  df_x_res_rcca2 <- tibble(ensembl_gene_id = colnames(X.mat),
                           weight = model_mat_rcca2$wX) %>% 
    left_join(df_modules_filt) %>% 
    mutate(abs_weight = abs(weight)) %>% 
    arrange(-abs_weight) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = unique(.$ensembl_gene_id)))
  
  df_x_res_rcca2 %>% 
    #filter(abs_weight > 0.0005) %>%
    ggplot(aes(x = ensembl_gene_id, y = abs_weight, color = module)) +
    geom_point() +
    scale_color_manual(values = c(df_x_res1 %>% 
                                    #head(100) %>% 
                                    dplyr::count(module, color) %>% 
                                    dplyr::select(-n) %>% 
                                    pull(color))
    ) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("absolute value of X coefficient weights")
  
  df_x_res_rcca2 %>% 
    slice_max(abs_weight, n = 30) %>% 
    arrange(weight) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = .$gene_symbol)) %>% 
    
    ggplot(aes(x = weight, y = gene_symbol, 
               fill = weight)) +
    geom_col() +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                         limits = c(-0.03, 0.03)) +
    ylab("") +
    theme_classic() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12)
    ) +
    ggtitle("Top 30 absolute gene weights")
  
  df_x_res_rcca2 %>% filter(module != 0) %>% 
    head(0.1*nrow(.)) %>% 
    count(module) %>% 
    left_join(df_modules_filt %>% 
                count(module) %>% 
                dplyr::rename("module_n" = "n")
    ) %>% 
    mutate(perc_mod = n/module_n * 100) %>% 
    arrange(-perc_mod)
  
  

# OVERLAP WITH NIRMALA RES
  load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res
  
  df_akula_res %>% 
    filter(id %in% c(
      
      df_x_res_rcca1 %>%
        #filter(abs_weight > 0.005) %>% 
        head(0.01*nrow(.)) %>%
        pull(ensembl_gene_id)
      
    ) & padj < 0.05
    
    ) %>% 
    
    left_join(df_modules_filt %>% dplyr::rename("id" = "ensembl_gene_id"))
    
  
  
# R PACKAGE BELOW ---------------------------------------------------------

  
# GRCCA functions -------------------------------------------------------------------
    
    
    ## FIX GRCCA FUNCTION
    GRCCA = function(X, Y, group1 = rep(1, ncol(X)), group2 = rep(1, ncol(Y)), lambda1 = 0, lambda2 = 0, mu1 = 0, mu2 = 0){
      n.comp = min(ncol(X), ncol(Y), nrow(X))
      
      #transform
      X.tr = GRCCA.tr(X, group1, lambda1, mu1)
      Y.tr = GRCCA.tr(Y, group2, lambda2, mu2)
      if(is.null(X.tr) | is.null(Y.tr)) return(invisible(NULL))
      
      #solve optimization problem
      sol = PRCCA(X.tr$mat, Y.tr$mat, X.tr$ind, Y.tr$ind, 1, 1)
      mod.cors = sol$mod.cors[1:n.comp]
      
      #inverse transform
      X.inv.tr = GRCCA.inv.tr(sol$x.coefs, group1, lambda1, mu1)
      x.coefs = as.matrix(X.inv.tr$coefs[,1:n.comp])
      rownames(x.coefs) = colnames(X)
      x.vars = sol$x.vars[,1:n.comp]
      rownames(x.vars) = rownames(X)
      Y.inv.tr = GRCCA.inv.tr(sol$y.coefs, group2, lambda2, mu2)
      y.coefs = as.matrix(Y.inv.tr$coefs[,1:n.comp])
      rownames(y.coefs) = colnames(Y)
      y.vars = sol$y.vars[,1:n.comp]
      rownames(y.vars) = rownames(Y)
      rho = diag(cor(x.vars,  y.vars))[1:n.comp]
      names(rho) = paste('can.comp', 1:n.comp, sep = '')
      return(list('n.comp' = n.comp, 'cors' = rho, 'mod.cors' = mod.cors, 'x.coefs' = x.coefs, 'x.vars' = x.vars, 'y.coefs' = y.coefs, 'y.vars' = y.vars))
    }
    
    # FROM SOURCE CODE
    
    GRCCA.tr = function(X, group, lambda, mu){
      if(is.null(group)){
        lambda = 0
        mu = 0
      } 
      if(lambda < 0){
        cat('please make', deparse(substitute(lambda)),'> 0\n')
        return(NULL)
      }
      if(mu < 0){
        cat('please make', deparse(substitute(mu)),'> 0\n')
        return(NULL)
      }
      index = NULL
      if(lambda > 0){
        #extend matrix
        group.names = unique(sort(group))
        ps = table(group)
        agg = aggregate(t(X), by = list(group), FUN = mean)
        X.mean = t(agg[, -1])
        colnames(X.mean) = agg[, 1] 
        if(mu == 0){
          mu = 1
          index = 1:ncol(X)
        } else {
          index = 1:(ncol(X) + ncol(X.mean))
        }
        X1 = 1/sqrt(lambda) * (X - X.mean[,group])
        X2 = scale(X.mean[,group.names], center = FALSE, scale = sqrt(mu/ps[group.names]))
        X = cbind(X1, X2)
      }
      return(list('mat' = X, 'ind' = index))
    }
    
    GRCCA.inv.tr = function(alpha, group, lambda, mu){
      if(is.null(group)){
        lambda = 0
        mu = 0
      } 
      if(lambda > 0){
        p = length(group)
        group.names = unique(sort(group))
        ps = table(group)
        alpha1 = alpha[1:p,]
        alpha2 = alpha[-(1:p),]
        agg = aggregate(alpha1, by = list(group), FUN = mean)
        alpha1.mean = agg[, -1]
        rownames(alpha1.mean) = agg[, 1]
        alpha1 = 1/sqrt(lambda) * (alpha1 - alpha1.mean[group,])
        if(mu == 0) mu = 1
        alpha2 = t(scale(t(alpha2[group.names,]), center = FALSE, scale = sqrt(mu*ps[group.names])))
        alpha = alpha1 + alpha2[group,]
      }
      return(list('coefs' = alpha))
    }
    
    
    
# hyperparameter optimization ---------------------------------------------


# FILTER X.MAT FOR NO GREY MODULE GENES
  X.mat <- X.mat[x_group != "M0"]
  x_group <- x_group[x_group != "M0"]
  
# SPLIT DATA
  set.seed(123)
  x_split <- initial_split(X.mat)
  x_training <- training(x_split)
  x_test <- testing(x_split)
  
  x_folds <- vfold_cv(x_training, v = 4)

# RUN GRCCA ON 5 VFOLDS TO GET 5 VECTORS OF X & Y COEFS FOR EACH SET OF PARAMETERS  
  df_hyper <- expand_grid(

    fold = 1:4,
    lambda1 = c(1 %o% 10^(-3:5)),
    mu1 = c(1 %o% 10^(-3:5))

  )

  doParallel::registerDoParallel()
  
  l_grcca_train <- pmap(
    
    .l = list(fold = df_hyper %>% pull(fold),
              lambda1 = df_hyper %>% pull(lambda1),
              mu1 = df_hyper %>% pull(mu1)
    ),
    
    .f = ~ list(
      
      print(paste0("fold = ", {..1}, ", lambda = ", {..2}, ", mu = ", {..3})),
      x_fold <- analysis(x_folds$splits[[ {..1} ]]),
      y_fold <- Y.mat[rownames(x_fold), ],
      GRCCA(X = x_fold,
            Y = y_fold,
            group1 = x_group,
            group2 = NULL,
            lambda1 = {..2}, # optimize
            lambda2 = 0,
            mu1 = {..3}, # optimize
            mu2 = 0) 
      
    )
    
  )
  
  save(l_grcca_train, file = "objects/20230228_grcca_train_noCovariates.RDS") # no grey genes?
  
  df_train_res <- map(l_grcca_train, 4) %>% 
    map(2) %>% 
    bind_rows %>% 
    bind_cols(df_hyper, .) %>% 
    pivot_longer(contains("comp"), names_to = "comp", values_to = "r")
  
  
# RUN ACROSS VALIDATION DATA  
  #load("objects/20230228_grcca_train.RDS")
  
  l_grcca_test <- pmap(
    
    .l = list(row = 1:nrow(df_hyper),
              fold = df_hyper %>% pull(fold),
              lambda1 = df_hyper %>% pull(lambda1),
              mu1 = df_hyper %>% pull(mu1)
    ),
    
    .f = ~ list(
      
      print(paste0("fold = ", {..2}, ", lambda = ", {..3}, ", mu = ", {..4})),
      x_val <- assessment(x_folds$splits[[ {..2} ]]),
      y_val <- Y.mat[rownames(x_val), ],
      
      cor(as.matrix(x_val) %*% l_grcca_train[[ {..1} ]][[4]]$x.coefs, 
          as.matrix(y_val) %*% l_grcca_train[[ {..1} ]][[4]]$y.coefs) %>% 
        
        as.data.frame %>% 
        rownames_to_column("comp") %>% 
        as_tibble() %>% 
        pivot_longer(contains("can.comp"), names_to = "comp2", values_to = "r") %>% 
        filter(comp == comp2) %>% 
        dplyr::select(-comp2) %>% 
        
        mutate(fold = {..2},
               lambda1 = {..3},
               mu1 = {..4}, 
               .before = 1)
      
    )
    
  )
  
  save(l_grcca_test, file = "objects/20230228_grcca_test_noCovariates.RDS")

  
  df_test_res <- map(l_grcca_test, 4) %>% 
    bind_rows  
  
  

# plot parameters & variate correlations ----------------------------------

  
# LOAD DATA
  
  # train res
  # load("objects/20230228_grcca_train.RDS") # l_grcca_train
  # 
  # df_train_res <- map(l_grcca_train, 4) %>%
  #   map(2) %>%
  #   bind_rows %>%
  #   bind_cols(df_hyper %>% filter(fold %in% 1:7), .) %>%
  #   pivot_longer(contains("comp"), names_to = "comp", values_to = "r")
  # 
  # rm(list = "l_grcca_train")
  # 
  # # test res
  # load("objects/20230228_grcca_test.RDS") # l_grcca_test
  # 
  # df_test_res <- map(l_grcca_test, 4) %>%
  #   bind_rows
  # 
  # rm(list = "l_grcca_test")



# COMBINE
  
  df_comp_cor <- df_train_res %>% 
    mutate(type = "train", .before = 1) %>% 
    bind_rows(df_test_res %>% 
                mutate(type = "test", .before = 1)) %>% 
    mutate(type = factor(type, levels = c("train", "test")))
  
  
# PLOT  
  df_comp_cor %>% 
    group_by(type, lambda1, mu1, comp) %>% 
    summarise(
      error = qnorm(0.975)*sd(r)/sqrt(5),
      mean = mean(r),
      low = mean - error,
      high = mean + error
    ) %>% 
    
    filter(comp == "can.comp1") %>% 
    
    ggplot(aes(x = log10(lambda1))) +
    geom_line(aes(y = mean, color = as.factor(mu1)), linewidth = 1.5) +
    geom_ribbon(aes(ymin = low, ymax = high, fill = as.factor(mu1)), alpha = 0.2) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_grid(type ~ comp) +
    ylab("correlation b/w variates") +
    ggtitle("training data correlation; lambda2 = 0 & mu2 = 0")
  
  df_comp_cor %>% 
    filter(type == "test") %>% 
    group_by(lambda1, mu1, comp) %>% 
    summarise(mean_r = mean(r)) %>% 
    pivot_wider(id_cols = c(lambda1, mu1), names_from = comp, values_from = mean_r) %>% 
    
    arrange(-can.comp1)
  
  
  
  
# run on test data ----------------------------------------------------

  
  set.seed(123)
  x_split <- initial_split(X.mat)
  x_training <- training(x_split)
  y_training <- Y.mat[rownames(x_training), ]
  
  x_test <- testing(x_split)
  y_test <- Y.mat[rownames(x_test), ]
  
# RUN GRCCA ON ENTIRE TRAINING SET
  grcca <- GRCCA(X = x_training,
                 Y = y_training,
                 group1 = x_group,
                 group2 = NULL,
                 lambda1 = 1, 
                 lambda2 = 0,
                 mu1 = 1000, #10
                 mu2 = 0)
  
# CHECK THAT IT WORKED
  as.matrix(x_test) %*% grcca$x.coefs %>% head
  grcca$x.vars %>% head
  
  as.matrix(y_test) %*% grcca$y.coefs %>% head
  grcca$y.vars %>% head
  
  cor(as.matrix(x_test) %*% grcca$x.coefs, as.matrix(y_test) %*% grcca$y.coefs)
  grcca$cors
  
# CALCULATE COEFFICIENTS ON TEST  
  cor(as.matrix(x_test) %*% grcca$x.coefs, as.matrix(y_test) %*% grcca$y.coefs)
  # 0.77!!
  
# explore results ---------------------------------------------------------

    
  
# EXPLORE COEFFICIENT VECTORS  
  df_y_coefs <- grcca$y.coefs %>% 
    as.data.frame %>% 
    rownames_to_column("dx") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    pivot_longer(starts_with("can"), names_to = "can_comp", values_to = "coef") %>% 
    mutate(can_comp = str_remove(can_comp, "can_comp")) %>% 
    mutate(dx = factor(dx, levels = c("bipolar", "mdd", "schizo", "age_death", "sex", "race", "anti_anxiety", "opioids")))
  
  df_x_coefs <- grcca$x.coefs %>% 
    as.data.frame %>% 
    rownames_to_column("group.gene") %>% 
    as_tibble() %>% 
    dplyr::rename("ensembl_gene_id" = "group.gene") %>% 
    pivot_longer(starts_with("can"), names_to = "can_comp", values_to = "coef") %>% 
    mutate(can_comp = str_remove(can_comp, "can.comp")) %>% 
    
    left_join(df_grcca_x %>% 
                dplyr::select(ensembl_gene_id, module),
              by = "ensembl_gene_id"
    ) %>% 
    distinct() %>% 
    group_by(can_comp, module) %>% 
    mutate(mean_coef = mean(coef)) %>% 
    ungroup 
  
  df_y_coefs %>% 
    filter(can_comp == 1) %>% 
    
    ggplot(aes(x = dx, y = coef, color = dx)) +
    geom_point(size = 3) +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    ggtitle("y coefficients of canonical component 1")
    facet_wrap( ~ can_comp, scales = "free")
  
# PLOT X COEFFICIENTS TOGETHER
  
    # assign color vector for modules
    col_vector <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"))
    col_vector <- col_vector[!col_vector %in% c("#666666")]
    
    df_x_coefs %>% 
      
      filter(module != 0) %>% 
      
      filter(can_comp == 1) %>% 
      arrange(module, coef) %>% 
      mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = .$ensembl_gene_id)) %>% 
      
      ggplot(aes(x = ensembl_gene_id, y = coef, color = module)) +
      geom_point(shape = 1) +
      geom_hline(yintercept = 0, lty = 2, color = "black") +
      #scale_fill_manual(values = col_vector) +
      scale_color_manual(values = col_vector) +
      theme(axis.text.x = element_blank()) +
      ggtitle("Gene weights for canonical component 1")
    
    
# Z-SCORE
    df_top_genes <- df_x_coefs %>% 
      filter(can_comp == 1) %>% 
      mutate(mean_coef = mean(coef),
             sd_coef = sd(coef),
             outlier = ifelse(abs(coef) > mean_coef + 3*sd_coef, 1, 0)) %>% 
      filter(outlier == 1) %>% 
      left_join(df_ensembl_to_symbol) %>% 
      mutate(gene_symbol = ifelse(is.na(gene_symbol), 
                                  ensembl_gene_id, 
                                  gene_symbol))
      
    
    mod_order <- df_top_genes %>% dplyr::count(module) %>% 
      left_join(
        
        df_modules_filt %>% 
          dplyr::count(module) %>% 
          dplyr::rename("module_n" = "n")
        
      ) %>% 
      mutate(perc_of_mod = n/module_n * 100) %>% 
      arrange(-perc_of_mod) %>% 
      
      pull(module) %>% unique

    # plot
    df_top_genes %>% 
      mutate(module = factor(module, levels = mod_order)) %>% 
      arrange(module, coef) %>% 
      mutate(gene_symbol = factor(gene_symbol, levels = unique(.$gene_symbol))) %>% 
      
      ggplot(aes(x = gene_symbol, y = coef, color = module)) +
      geom_point(size = 2) +
      scale_color_manual(values = col_vector) +
      facet_wrap(vars(module), scales = "free") +
      ggtitle("Genes with abs(coef) > mean + 3*SD") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
            strip.text = element_text(size = 10),
            plot.title = element_text(size = 10),
            legend.position = "none",
            axis.title.x = element_blank())
      
     
    df_top_genes %>% 
      filter(module == 13) %>% 
      arrange(-abs(coef)) %>% 
      left_join(df_ensembl_to_symbol)
  
  
    
# IDENTIFY GENES OF INTEREST
  
  grcca_genes <- df_x_coefs %>% 
    filter(module %in% c(4, 9, 14, 16) & can_comp == 1 &
             (coef > upper | coef < lower)) %>% 
    pull(ensembl_gene_id)
  
  
  # get gene info from ensembl database
  hsapiens_ensembl <- useEnsembl(biomart = "genes",
                                 dataset = "hsapiens_gene_ensembl")
  
  listAttributes(hsapiens_ensembl) %>% 
    as_tibble %>% 
    filter(grepl("def", name))
    
  m_grcca_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "definition_1006", "gene_biotype"),
        filters = "ensembl_gene_id",
        values = grcca_genes,
        mart = hsapiens_ensembl)
  
  m_grcca_genes %>% 
    as_tibble() %>% 
    filter(definition_1006 != "") %>% 
    dplyr::select(-definition_1006) %>% 
    distinct() %>% 
    left_join(df_x_coefs %>% 
                filter(module %in% c(4, 9, 14, 16) & can_comp == 1 &
                         (coef > upper | coef < lower))) %>% 
    dplyr::select(module, external_gene_name) %>% 
    arrange(module)
  
# EXPLORE LVs
  
  df_latent_vars <- grcca$x.vars %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    pivot_longer(starts_with("can"), names_to = "can_comp", values_to = "coef") %>% 
    mutate(can_comp = str_remove(can_comp, "can_comp")) %>%
    mutate(dim = "x", .before = 1) %>% 
    
    bind_rows(
      grcca$y.vars %>% 
        as.data.frame %>% 
        rownames_to_column("sample") %>% 
        as_tibble() %>% 
        clean_names() %>% 
        pivot_longer(starts_with("can"), names_to = "can_comp", values_to = "coef") %>% 
        mutate(can_comp = str_remove(can_comp, "can_comp")) %>% 
        mutate(dim = "y", .before = 1)
    )
  
  df_latent_vars %>% 
    pivot_wider(id_cols = c(sample, can_comp), names_from = dim, values_from = coef) %>% 
    
    ggplot(aes(x = x, y = y)) +
    geom_point() +
    facet_wrap( ~ can_comp, scales = "free")
  
