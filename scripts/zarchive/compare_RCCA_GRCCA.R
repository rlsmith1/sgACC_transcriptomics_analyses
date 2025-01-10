
# libraries ---------------------------------------------------------------


  library(tidyverse)
  library(patchwork)
  library(janitor)
  library(RColorBrewer)
  library(RCCA)
  library(biomaRt)
  library(tidymodels)  
  library(ggrepel)
  library(raveio)





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
  

# GRCCA: lambda = 0.9999; mu = 0.2 ----------------------------------------


# RES 1
  cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_19kgenes_covarOpioidsAnx_regressBrainweightAgedeath/")
  model_mat1 <- read_mat(paste0(cca_dir, "framework/grcca_holdout1-0.20/res/level1/model_1.mat"))
  
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
    stat_cor(label.y = 0.00, label.x = -4.5) +
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
    left_join(df_ensembl_to_symbol) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    arrange(module) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = .$ensembl_gene_id))
  
  # plot all weights
  df_col <- df_modules_filt %>% count(module, color)
  colors <- df_col %>% pull(color)
  names(colors) <- df_col %>% pull(module)

  df_x_res1 %>% 
    ggplot(aes(x = ensembl_gene_id, y = weight)) +
    geom_col(aes(fill = module)) +
    scale_fill_manual(values = colors) +
    theme(axis.text.x = element_blank())
  
  # plot top-weighted genes
  df_x_res1 %>% 
    slice_max(abs_weight, n = 30) %>% 
    arrange(weight) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = .$gene_symbol)) %>% 
    ggplot(aes(x = weight, y = gene_symbol, fill = weight)) +
    geom_col() +
    #coord_flip() +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                         limits = c(-0.02, 0.02)) +
    ylab("") +
    theme_classic()
    
  # top modules?
  df_x_res1 %>% 
    group_by(module) %>% 
    summarise(mean_weight = mean(abs_weight)) %>% 
    arrange(-mean_weight) %>% 
    print(n = nrow(.))
  
  
  
# RES 2
  cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_19kgenes_covarOpioidsAnx_regressBrainweightAgedeath/")
  model_mat2 <- read_mat(paste0(cca_dir, "framework/grcca_holdout1-0.20/res/level2/model_1.mat"))
  
# LEVEL 2 X-Y COR
  df_lvs2 <- tibble(lvx = as.matrix(X.mat) %*% model_mat2$wX,
                    lvy = as.matrix(Y.mat) %*% model_mat2$wY)
  
  df_lvs2 %>% 
    ggplot(aes(x = lvx, y = lvy)) +
    geom_point() +
    geom_smooth(method = lm, lty = 2) +
    xlab("X latent variable (x matrix • x weights)") +
    ylab("Y latent variable (y matrix • y weights)") +
    #stat_regline_equation() +
    stat_cor(label.y = 0.07, label.x = -3) +
    ggtitle("Gene x-y latent variable correlation")
  
  
# LEVEL 2 Y   
  df_y_res2 <- tibble(covariate = colnames(Y.mat), 
                      weight = model_mat2$wY) %>% 
    mutate(covariate = factor(covariate, levels = c("bipolar", "mdd", "schizo", 
                                                    "opioids", "anti_anxiety")))
  
  df_y_res2 %>% 
    ggplot(aes(x = covariate, y = weight)) +
    geom_point(size = 5, aes(color = weight)) +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    scale_color_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                          limits = c(-0.1, 0.1)) +
    ylim(c(-0.1, 0.1)) +
    xlab("") +
    ggtitle("Y coefficient weights") +
    theme(legend.position = "none") # opioid effect
  
# LEVEL 2 X  
  df_x_res2 <- tibble(ensembl_gene_id = colnames(X.mat),
                      weight = model_mat2$wX) %>% 
    left_join(df_modules_filt) %>% 
    mutate(abs_weight = abs(weight)) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    arrange(module) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = .$ensembl_gene_id))
  
  # plot all weights
  df_col <- df_modules_filt %>% count(module, color)
  colors <- df_col %>% pull(color)
  names(colors) <- df_col %>% pull(module)
  
  df_x_res2 %>% 
    ggplot(aes(x = ensembl_gene_id, y = weight)) +
    geom_col(aes(fill = module)) +
    scale_fill_manual(values = colors) +
    theme(axis.text.x = element_blank())
  
  # plot top-weighted genes
  df_x_res2 %>% 
    slice_max(abs_weight, n = 30) %>% 
    arrange(weight) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = .$gene_symbol)) %>% 
    ggplot(aes(x = weight, y = gene_symbol, fill = weight)) +
    geom_col() +
    #coord_flip() +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                         limits = c(-0.02, 0.02)) +
    ylab("") +
    theme_classic()
  
  # top modules?
  df_x_res2 %>% 
    group_by(module) %>% 
    summarise(mean_weight = mean(abs_weight)) %>% 
    arrange(-mean_weight) %>% 
    print(n = nrow(.))
  

# RES 3
  model_mat3 <- read_mat(paste0(cca_dir, "framework/grcca_holdout1-0.20/res/level3/model_1.mat"))
  
  # LEVEL 2 X-Y COR
  df_lvs3 <- tibble(lvx = as.matrix(X.mat) %*% model_mat3$wX,
                    lvy = as.matrix(Y.mat) %*% model_mat3$wY)
  
  df_lvs3 %>% 
    ggplot(aes(x = lvx, y = lvy)) +
    geom_point() +
    geom_smooth(method = lm, lty = 2) +
    xlab("X latent variable (x matrix • x weights)") +
    ylab("Y latent variable (y matrix • y weights)") +
    #stat_regline_equation() +
    stat_cor(label.y = 0.07, label.x = -3) +
    ggtitle("Gene x-y latent variable correlation")
  
  
  # LEVEL 3 Y   
  df_y_res3 <- tibble(covariate = colnames(Y.mat), 
                      weight = model_mat3$wY) %>% 
    mutate(covariate = factor(covariate, levels = c("bipolar", "mdd", "schizo", 
                                                    "opioids", "anti_anxiety")))
  
  df_y_res3 %>% 
    ggplot(aes(x = covariate, y = weight)) +
    geom_point(size = 5, aes(color = weight)) +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    scale_color_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                          limits = c(-0.1, 0.1)) +
    ylim(c(-0.1, 0.1)) +
    xlab("") +
    ggtitle("Y coefficient weights") +
    theme(legend.position = "none") # opioid effect
  
  # LEVEL 3 X  
  df_x_res3 <- tibble(ensembl_gene_id = colnames(X.mat),
                      weight = model_mat3$wX) %>% 
    left_join(df_modules_filt) %>% 
    mutate(abs_weight = abs(weight)) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    arrange(module) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = .$ensembl_gene_id))
  
  # plot all weights
  df_x_res3 %>% 
    ggplot(aes(x = ensembl_gene_id, y = weight)) +
    geom_col(aes(fill = module)) +
    scale_fill_manual(values = colors) +
    theme(axis.text.x = element_blank())
  
  # plot top-weighted genes
  df_x_res3 %>% 
    slice_max(abs_weight, n = 30) %>% 
    arrange(weight) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = .$gene_symbol)) %>% 
    ggplot(aes(x = weight, y = gene_symbol, fill = weight)) +
    geom_col() +
    #coord_flip() +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                         limits = c(-0.023, 0.02)) +
    ylab("") +
    theme_classic()
  
  # top modules?
  df_x_res3 %>% 
    group_by(module) %>% 
    summarise(mean_weight = mean(abs_weight)) %>% 
    arrange(-mean_weight) %>% 
    print(n = nrow(.))
  
  
  

# SCZ genome overlap ------------------------------------------------------

  
# SCZ DRUGGABLE GENOME  
  
  df_scz_genome <- readxl::read_excel("data/eva/lago_bahn_2022_scz_genome.xlsx") %>% 
    row_to_names(1) %>% 
    clean_names
  
  scz_genes <- df_scz_genome$risk_gene_hgnc %>% unique()
  
  gene_universe <- df_x_res1$gene_symbol %>% unique
  
  scz_genes_in_universe <- intersect(gene_universe, scz_genes) 
  
  top_grcca_genes <- df_x_res1 %>%
    arrange(-abs_weight) %>% 
    head(20) %>% 
    #head(0.001*nrow(.)) %>% 
    pull(gene_symbol) # significant
  
  overlap_genes <- intersect(top_grcca_genes, scz_genes_in_universe)
  
  phyper(q = length(overlap_genes) - 1, 
         m = length(scz_genes_in_universe), 
         n = length(gene_universe) - length(scz_genes_in_universe), 
         k = length(top_grcca_genes), 
         lower.tail = FALSE, log.p = FALSE)
  

# SCZ GWAS  
  library(gprofiler2)
  df_scz_gwas_genes <- readxl::read_excel("data/eva/trubetskoy_et_al_2022_scz_GWAS.xls") %>% 
    clean_names %>%
    pull(top_index) %>% 
    gsnpense %>% 
    as_tibble()
  
  scz_gwas_genes <- df_scz_gwas_genes %>% pull(gene_names) %>% unique
  gene_universe <- df_x_res1 %>% pull(gene_symbol) %>% unique
  scz_gwas_genes_in_universe <- intersect(gene_universe, scz_gwas_genes) 
  
  top_grcca_genes <- df_x_res1 %>%
    arrange(-abs_weight) %>% 
    head(10) %>% 
    pull(gene_symbol)
  
  overlap_genes <- intersect(top_grcca_genes, scz_gwas_genes_in_universe)
  
  q <- length(overlap_genes) - 1
  m <- length(trubetskoy_genes_in_universe)
  n <- length(gene_universe) - m
  k <- length(top_cca_genes)
  phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  
# OVERLAP WITH NIRMALA RES
  load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res
  
  df_akula_res %>% 
    filter(id %in% c(
      
      df_x_res1 %>%
        head(200) %>%
        pull(ensembl_gene_id)
      
    ) & padj < 0.05
    
    ) %>% 
    
    left_join(df_modules_filt %>% dplyr::rename("id" = "ensembl_gene_id"))
  
