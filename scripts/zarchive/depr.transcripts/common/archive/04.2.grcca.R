
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





# set theme for plots -----------------------------------------------------


  theme_set(theme_bw() +
              theme(plot.title = element_text(size = 18),
                    axis.title = element_text(size = 15),
                    axis.text = element_text(size = 12),
                    strip.text = element_text(size = 15),
                    legend.title = element_text(size = 15),
                    legend.text = element_text(size = 12)))
  



# data --------------------------------------------------------------------


base_dir <- "/Users/work/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/projects/WGCNA/acsg_wgcna/"
prefix <- "20230314_185samples_26kCOMMONtranscripts_vst_qSVA123567_MP_RNAbatch_Race_resids"

# LOAD OBJECTS  
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress_common
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT3_MODS.RDS")) # df_modules
load(paste0(base_dir, "objects/transcript_gene_go_term.RDS")) # df_transcript_gene_go_term

# IDENTIFY MODULE SET OF INTEREST  
df_modules_filt <- df_modules %>% 
  filter(cut_height == 0.99 & min_size == 50) %>% 
  unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
  arrange(mod_set, module)

# COVARIATES
load("data/covariates/185_all_covariates_clean.Rdata")

# GENE SYMBOL TO ENSEMBL GENE ID MAPPING
load("objects/transcript_gene_go_term.RDS") # df_transcript_gene_go_term
load("objects/df_hsapiens_genome.RDS")
  
  
# define matrices ---------------------------------------------------------------
  
  
## SELECT COVARIATES BASED ON 4.1 RES: 
# anticholinergics, opioids, anti-epileptics, age_death


# DEFINE MATRICES  
  df_grcca_x <- df_vsd_regress_common %>% 
    dplyr::select(ensembl_transcript_id, sample, resids) %>% 
    left_join(df_modules_filt) %>% 
    arrange(module) 
  
  X.mat <- df_grcca_x %>% 
    pivot_wider(id_cols = sample, names_from = ensembl_transcript_id, values_from = resids) %>% 
    as.data.frame %>% 
    column_to_rownames("sample")
  
  x_group <- df_grcca_x %>% 
    dplyr::select(ensembl_transcript_id, module) %>% 
    distinct() %>% 
    pull(module) %>% 
    paste0("M", .)
  
  # filter x for no grey module genes
  #X.mat <- X.mat[x_group != "M0"]
  #x_group <- x_group[x_group != "M0"]
  
  # remove M for Agoston's toolbox
  x_group <- str_remove(x_group, "M") %>% as.numeric
  
  YC.mat <- df_covariates_all_clean %>% 
    dplyr::select(sample, opioids, anticholinergics, anti_epileptics, age_death) %>%  
    mutate(control = ifelse(grepl("control", sample), 1, 0),
           bipolar = ifelse(grepl("bipolar", sample), 1, 0),
           mdd = ifelse(grepl("mdd", sample), 1, 0),
           schizo = ifelse(grepl("schizo", sample), 1, 0),
           
           anti_epileptics = ifelse(anti_epileptics == "positive", 1, 0),
           anticholinergics = ifelse(anticholinergics == "positive", 1, 0),
           opioids = ifelse(opioids == "positive", 1, 0)#,
    ) %>% 
    dplyr::select(sample, control, bipolar, mdd, schizo, 
                  opioids, anti_epileptics, anticholinergics, age_death) %>%
    column_to_rownames("sample") %>% 
    as.data.frame
  
  Y.mat <- YC.mat %>% dplyr::select(-c(control, age_death)) # remove control column so columns aren't linearly equivalent
  
  C.mat <- YC.mat %>% dplyr::select(-c(control, bipolar, mdd, schizo, opioids, anticholinergics, anti_epileptics))
  
  # EXPORT FOR AGOSTON'S TOOLBOX
  # cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_26kCOMMONtranscripts_covarOpioidsCholEpilep_regressAgedeath/")
  # dir.create(cca_dir)
  # cca_data_dir <- paste0(cca_dir, "data/")
  # dir.create(cca_data_dir)
  # 
  # write.table(X.mat, file = paste0(cca_data_dir, "X.txt"))
  # write.table(Y.mat, file = paste0(cca_data_dir, "Y.txt"))
  # write.table(C.mat, file = paste0(cca_data_dir, "C.txt"))
  # write.table(x_group, file = paste0(cca_data_dir, "XGroup.txt"))
  
  cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_26kCOMMONtranscripts_covarOpioidsCholEpilep_regressAgedeath/")
  
  df_labels_x <- tibble(ensembl_transcript_id = colnames(X.mat)) %>% 
    left_join(df_modules_filt) %>% 
    left_join(df_hsapiens_genome %>% 
                dplyr::select(transcript_id, transcript_name) %>% 
                distinct %>% 
                dplyr::rename("ensembl_transcript_id" = "transcript_id")) %>% 
    mutate(Label = row_number(), .before = 1) %>% 
    dplyr::select(Label, module, ensembl_transcript_id, transcript_name) %>% 
    dplyr::rename_all(~c("Label", "Category", "Ensembl_id", "HGNC_id")) %>% 
    mutate(HGNC_id = ifelse(is.na(HGNC_id), Ensembl_id, HGNC_id))
  
  df_labels_y <- tibble(Label = colnames(Y.mat))
  
  write.csv(df_labels_x, paste0(cca_dir, "data/LabelsX.csv"), row.names = FALSE)
  write.csv(df_labels_y, paste0(cca_dir, "data/LabelsY.csv"), row.names = FALSE)
  

  
# GRCCA results -----------------------------------------------------------
  
  
### LEVEL 1
  cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_26kCOMMONtranscripts_covarOpioidsCholEpilep_regressAgedeath/")
  model_mat1 <- read_mat(paste0(cca_dir, "framework/grcca_holdout1-0.20_subsamp5-0.20/res/level1/model_1.mat"))
  param_mat1 <- read_mat(paste0(cca_dir, "framework/grcca_holdout1-0.20_subsamp5-0.20/res/level1/param_1.mat"))
  lambda <- param_mat1$`param/L2x`
  mu <- param_mat1$`param/groupL2x`
  
# LEVEL 1 X-Y COR
  df_lvs1 <- tibble(lvx = as.matrix(X.mat) %*% model_mat1$wX,
                    lvy = as.matrix(Y.mat) %*% model_mat1$wY)
  
  df_lvs1 %>% 
    ggplot(aes(x = lvx, y = lvy)) +
    geom_point() +
    geom_smooth(method = lm, lty = 2, color = "midnightblue") +
    xlab("X latent variable (x matrix • x weights)") +
    ylab("Y latent variable (y matrix • y weights)") +
    #stat_regline_equation() +
    stat_cor(label.x = -.5, label.y = 0.1) +
    ggtitle("Rare transcript x-y latent variable correlation")
  
# LEVEL 1 Y   
  df_y_res1 <- tibble(covariate = colnames(Y.mat), 
                      weight = model_mat1$wY) %>% 
    mutate(covariate = factor(covariate, levels = c("bipolar", "mdd", "schizo", 
                                                    "opioids", "anticholinergics", "anti_epileptics")))
  
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
  df_x_res1 <- tibble(transcript_id = colnames(X.mat),
                      weight = model_mat1$wX) %>% 
    left_join(df_hsapiens_genome) %>% 
    left_join(df_modules_filt %>% 
                dplyr::rename("transcript_id" = "ensembl_transcript_id")) %>% 
    mutate(abs_weight = abs(weight)) %>% 
    arrange(-abs_weight) %>% 
    mutate(transcript_id = factor(transcript_id, levels = unique(.$transcript_id))) 
  
  df_x_res1 %>% 
    #filter(abs_weight > 0.0005) %>%
    ggplot(aes(x = transcript_id, y = abs_weight, color = module)) +
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
    mutate(transcript_name = ifelse(is.na(transcript_name), 
                                    as.character(transcript_id), 
                                    transcript_name)) %>% 
    
    dplyr::select(transcript_id, transcript_name, weight, abs_weight) %>% distinct %>% 
    slice_max(abs_weight, n = 20) %>% 
    arrange(weight) %>% 
    mutate(transcript_name = factor(transcript_name, levels = .$transcript_name)) %>% 
    
    ggplot(aes(x = weight, y = transcript_name, fill = weight)) +
    geom_col() +
    #coord_flip() +
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
    ggtitle("Top 20 absolute transcript weights")
  
  
### LEVEL 2
  cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_26kCOMMONtranscripts_covarOpioidsCholEpilep_regressAgedeath/")
  model_mat2 <- read_mat(paste0(cca_dir, "framework/grcca_holdout1-0.20_subsamp5-0.20/res/level2/model_1.mat"))
  param_mat2 <- read_mat(paste0(cca_dir, "framework/grcca_holdout1-0.20_subsamp5-0.20/res/level2/param_1.mat"))
  lambda <- param_mat2$`param/L2x`
  mu <- param_mat2$`param/groupL2x`
  
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
    stat_cor(label.x = -15, label.y = 0.07) +
    ggtitle("Rare transcript x-y latent variable correlation")
  
# LEVEL 2 Y   
  df_y_res2 <- tibble(covariate = colnames(Y.mat), 
                      weight = model_mat2$wY) %>% 
    mutate(covariate = factor(covariate, levels = c("bipolar", "mdd", "schizo", 
                                                    "opioids", "anti_epileptics", "anticholinergics")))
  
  df_y_res2 %>% 
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
  df_x_res2 <- tibble(transcript_id = colnames(X.mat),
                      weight = model_mat2$wX) %>% 
    left_join(df_hsapiens_genome) %>% 
    left_join(df_modules_filt %>% 
                dplyr::rename("transcript_id" = "ensembl_transcript_id")) %>% 
    mutate(abs_weight = abs(weight)) %>% 
    arrange(-abs_weight) %>% 
    mutate(transcript_id = factor(transcript_id, levels = unique(.$transcript_id))) 
  
  df_x_res2 %>% 
    #filter(abs_weight > 0.0005) %>%
    ggplot(aes(x = transcript_id, y = abs_weight, color = module)) +
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
  
  df_x_res2 %>% 
    mutate(transcript_name = ifelse(is.na(transcript_name), 
                                    as.character(transcript_id), 
                                    transcript_name)) %>% 
    
    dplyr::select(transcript_id, transcript_name, weight, abs_weight) %>% distinct %>% 
    slice_max(abs_weight, n = 20) %>% 
    arrange(weight) %>% 
    mutate(transcript_name = factor(transcript_name, levels = .$transcript_name)) %>% 
    
    ggplot(aes(x = weight, y = transcript_name, fill = weight)) +
    geom_col() +
    #coord_flip() +
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
    ggtitle("Top 20 absolute transcript weights")
  
  
  
# RCCA results ------------------------------------------------------------
  
  
################## LEVEL 1    
  
# LEVEL 1 RES    
  cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_26kCOMMONtranscripts_covarOpioidsCholEpilep_regressAgedeath/")
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
    stat_cor(label.y = 0.02, label.x = -30) +
    ggtitle("Transcript x-y latent variable correlation")
  
# LEVEL 1 Y   
  df_y_res_rcca1 <- tibble(covariate = colnames(Y.mat), 
                           weight = model_mat_rcca1$wY) %>% 
    mutate(covariate = factor(covariate, levels = c("bipolar", "mdd", "schizo", 
                                                    "anti_epileptics", "opioids", "anticholinergics")))
  
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
  df_x_res_rcca1 <- tibble(ensembl_transcript_id = colnames(X.mat),
                           weight = model_mat_rcca1$wX) %>% 
    left_join(df_modules_filt) %>% 
    left_join(df_hsapiens_genome %>% 
                dplyr::select(transcript_id, transcript_name) %>% 
                dplyr::rename("ensembl_transcript_id" = "transcript_id") %>% 
                distinct()) %>% 
    
    mutate(abs_weight = abs(weight)) %>% 
    arrange(-abs_weight) %>% 
    mutate(transcript_name = ifelse(is.na(transcript_name), ensembl_transcript_id, transcript_name)) %>% 
    mutate(transcript_name = factor(transcript_name, levels = unique(.$transcript_name)))
  
  df_x_res_rcca1 %>% 
    ggplot(aes(x = transcript_name, y = abs_weight, color = module)) +
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
  
  # plot top n
  df_col <- df_modules_filt %>% 
    dplyr::select(module, color) %>% 
    distinct()
  colors <- df_col %>% pull(color)
  names(colors) <- df_col %>% pull(module)
  
  df_x_res_rcca1 %>% 
    slice_max(abs_weight, n = 30) %>% 
    arrange(weight) %>% 
    mutate(transcript_name = factor(transcript_name, levels = unique(.$transcript_name))) %>% 
  
    ggplot(aes(x = weight, y = transcript_name, 
               fill = weight)) +
    geom_col() +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                         limits = c(-0.021, 0.022)) +
    ylab("") +
    theme_classic() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12)
    ) +
    ggtitle("Top 30 absolute transcript weights")
  
  df_x_res_rcca1 %>% filter(module != 0) %>% 
    head(0.1*nrow(.)) %>% 
    count(module) %>% 
    left_join(df_modules_filt %>% 
                count(module) %>% 
                dplyr::rename("module_n" = "n")
    ) %>% 
    mutate(perc_mod = n/module_n * 100) %>% 
    arrange(-perc_mod)
  
  

# contextualize results ---------------------------------------------------


# VARIANT TYPE
  
  # top RCCA transcripts
  df_rcca_transcript_biotype <- df_hsapiens_genome %>% 
    dplyr::filter(transcript_id %in% c(df_x_res_rcca1 %>% 
                                         head(0.01*nrow(.)) %>% 
                                         pull(ensembl_transcript_id))
    ) %>% 
    dplyr::select(transcript_id, transcript_biotype) %>% 
    distinct() %>% 
    count(transcript_biotype) %>% 
    arrange(-n) %>% 
    mutate(transcript_biotype = factor(transcript_biotype, 
                                       levels = .$transcript_biotype))
  
  p_rcca_biotype <- df_rcca_transcript_biotype %>% 
    ggplot(aes(x = transcript_biotype, y = n)) +
    geom_col(fill = "seagreen", alpha = 0.7) +
    geom_text(aes(label = n), vjust = -0.2) +
    labs(x = "") +
    ggtitle("Biotype of top 1% weighted transcript (n = 262)") +
    theme(axis.text.x = element_text(angle = 27, hjust = 1))
  
  # rare transcripts
  df_common_transcript_biotype <- df_hsapiens_genome %>% 
    dplyr::filter(transcript_id %in% c(df_x_res_rcca1 %>% 
                                         pull(ensembl_transcript_id))
    ) %>% 
    dplyr::select(transcript_id, transcript_biotype) %>% 
    distinct() %>% 
    filter(transcript_biotype %in% df_rcca_transcript_biotype$transcript_biotype) %>% 
    count(transcript_biotype) %>% 
    mutate(transcript_biotype = factor(transcript_biotype, 
                                       levels = df_rcca_transcript_biotype$transcript_biotype))
  
  p_rare_biotype <- df_common_transcript_biotype %>% 
    ggplot(aes(x = transcript_biotype, y = n)) +
    geom_col(fill = "midnightblue", alpha = 0.7) +
    geom_text(aes(label = n), vjust = -0.2) +
    labs(x = "") +
    ggtitle("Common transcript biotype (n = 26237)") +
    theme(axis.text.x = element_text(angle = 27, hjust = 1))
  
  p_rcca_biotype + p_rare_biotype
  
  # hypergeometric overlap
  df_hypergeometric_rare_tx <- df_common_transcript_biotype %>% 
    dplyr::rename("total_in_genome" = "n") %>% 
    left_join(
      
      df_rcca_transcript_biotype %>% 
        dplyr::rename("n_in_top_cca_transcripts" = "n")
      
    ) %>% 
    filter(!is.na(n_in_top_cca_transcripts)) %>% 
    
    mutate(
      
      q = n_in_top_cca_transcripts - 1,
      m = total_in_genome,
      n = nrow(df_x_res_rcca1) - m,
      k = 0.01*nrow(df_x_res_rcca1),
      p_value = phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE),
      p_adj = p.adjust(p_value, method = "BH")
      
    )
  
  df_hypergeometric_rare_tx %>% filter(p_value < 0.05)
  
# DRUGGABLE SCZ GENOME OVERLAP
  
  # data
  load("objects/df_hsapiens_genome.RDS") # df_hsapiens_genome
  df_scz_genome <- readxl::read_excel("data/eva/lago_bahn_2022_scz_genome.xlsx") %>% 
    row_to_names(1) %>% 
    clean_names
  
  # hypergeometric test
  lago_bahn_risk_genes <- df_scz_genome$risk_gene_hgnc %>% unique
  
  gene_universe <- df_hsapiens_genome %>% 
    filter(transcript_id %in% colnames(X.mat)) %>% 
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>% 
    pull(gene_name) %>% unique # 13506
  
  lago_bahn_risk_genes_in_universe <- intersect(gene_universe, lago_bahn_risk_genes) # 567
  
  top_genes_rcca1 <- df_x_res_rcca1 %>%
    head(0.01*nrow(.)) %>%
    left_join(df_hsapiens_genome %>% 
                dplyr::select(transcript_id, gene_name, gene_id) %>% 
                distinct() %>% 
                dplyr::rename("ensembl_transcript_id" = "transcript_id")) %>% 
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>% 
    pull(gene_name) # 262
  
  overlap_genes1 <- intersect(top_genes_rcca1, lago_bahn_risk_genes_in_universe) # 10
  
  q <- length(overlap_genes1) - 1
  m <- length(lago_bahn_risk_genes_in_universe)
  n <- length(gene_universe) - m
  k <- length(top_genes_rcca1)
  
  phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  
  # plot overlap genes
  df_x_res1 %>% 
    filter(gene_symbol %in% overlap_genes & abs_weight > 0.005) %>% 
    arrange(weight) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = unique(.$gene_symbol))) %>% 
    ggplot(aes(x = weight, y = gene_symbol, fill = weight)) +
    geom_col() +
    #coord_flip() +
    scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                         limits = c(-0.0151, 0.0151)) +
    ylab("") +
    theme_classic() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 15),
          legend.position = "none"
    ) +
    ggtitle("SCZ druggable targets; top 1% genes")
  
  # TRANSDIAGNOSTIC
  df_pgc_gwas <- readxl::read_excel("data/PGC_transdiagnostic_GWAS.xlsx") %>% 
    row_to_names(2) %>% 
    clean_names
  
  # pgc_genes <- df_pgc_gwas %>% 
  #   dplyr::select(candidate_gene) %>% 
  #   mutate(candidate_gene = gsub("\\s*\\([^\\)]+\\)", "", candidate_gene)) %>% 
  #   filter(candidate_gene != "-") %>% 
  #   separate(candidate_gene, into = paste0('gene', seq_len(max(str_count(.$candidate_gene, ';')+1))), 
  #            sep = ";", fill = "right") %>% 
  #   pivot_longer(1:ncol(.), names_to = "gene", values_to = "gene_symbol") %>% 
  #   filter(!is.na(gene_symbol)) %>% 
  #   pull(gene_symbol) %>% 
  #   unique()
  
  pgc_genes <- df_pgc_gwas %>% 
    dplyr::select(nearest_gene_100kb) %>% 
    mutate(nearest_gene_100kb = gsub("\\s*\\([^\\)]+\\)", "", nearest_gene_100kb)) %>% 
    filter(nearest_gene_100kb != "-") %>% 
    separate(nearest_gene_100kb, 
             into = paste0('gene', seq_len(max(str_count(.$nearest_gene_100kb, ',')+1))), 
             sep = ",", fill = "right") %>% 
    pivot_longer(1:ncol(.), names_to = "gene", values_to = "gene_symbol") %>% 
    filter(!is.na(gene_symbol)) %>% 
    pull(gene_symbol) %>% 
    unique()
  
  pgc_genes_in_universe <- intersect(gene_universe, pgc_genes) 
  
  top_genes_rcca1 <- df_x_res_rcca1 %>%
    head(0.01*nrow(.)) %>%
    left_join(df_hsapiens_genome %>% 
                dplyr::select(transcript_id, gene_name, gene_id) %>% 
                distinct() %>% 
                dplyr::rename("ensembl_transcript_id" = "transcript_id")) %>% 
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>% 
    pull(gene_name) # 262
  
  overlap_genes1 <- intersect(top_genes_rcca1, pgc_genes_in_universe) # 10
  
  q <- length(overlap_genes1) - 1
  m <- length(pgc_genes_in_universe)
  n <- length(gene_universe) - m
  k <- length(top_genes_rcca1)
  phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  
  
  # PLOT CHROMOSOME REGION
  
  # get gene info from ensembl database
  hsapiens_ensembl <- useEnsembl(biomart = "genes",
                                 dataset = "hsapiens_gene_ensembl")
  
  listAttributes(hsapiens_ensembl) %>% 
    tibble %>% 
    filter(grepl("chrom", name))
  
  df_chromosomes <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                          filters = "hgnc_symbol",
                          values = overlap_genes1,
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
                dplyr::rename("transcript_id" = "ensembl_transcript_id") %>% 
                left_join(df_hsapiens_genome %>% 
                            dplyr::select(transcript_id, gene_id, gene_name) %>% 
                            distinct() %>% 
                            filter(gene_name %in% overlap_genes1)) %>% 
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
  
  
# OVERLAP WITH NIRMALA RES
  load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res
  
  df_akula_res %>% 
    filter(id %in% c(
      
      df_x_res_rcca1 %>%
        head(0.01*nrow(.)) %>%
        pull(ensembl_transcript_id)
      
    ) & padj < 0.05
    
    ) %>% 
    
    left_join(df_modules_filt %>% dplyr::rename("id" = "ensembl_transcript_id"))
  
  
  
  
# PCA-CCA results -----------------------------------------------
  
  
  
  #cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_26kCOMMONtranscripts_regressOpioidsSuicideBrainweightAgedeath/")
  #cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_26kCOMMONtranscripts_covarOpioidsSuicide_regressBrainweightAgedeath/")
  #cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_26kCOMMONtranscripts_regressOpioidsCholEpilepAgedeath/")
  cca_dir <- paste0(base_dir, "RCCA_toolkit/s185_26kCOMMONtranscripts_covarOpioidsChol_regressEpilepAgedeath/")
  model_mat <- read_mat(paste0(cca_dir, "framework/cca_holdout1-0.20_subsamp5-0.20/res/level1/model_1.mat"))
  
# LEVEL 1 X-Y COR
  df_lvs <- tibble(lvx = as.matrix(X.mat) %*% model_mat$wX,
                   lvy = as.matrix(Y.mat) %*% model_mat$wY)
  
  df_lvs %>% 
    ggplot(aes(x = lvx, y = lvy)) +
    geom_point() +
    geom_smooth(method = lm, lty = 2) +
    xlab("X latent variable (x matrix • x weights)") +
    ylab("Y latent variable (y matrix • y weights)") +
    #stat_regline_equation() +
    stat_cor(label.y = 0.025, label.x = -0.3) +
    ggtitle("Common transcript x-y latent variable correlation")
  
  
# LEVEL 1 Y   
  df_y_res <- tibble(covariate = colnames(Y.mat), 
                     weight = model_mat$wY) %>% 
    mutate(covariate = factor(covariate, levels = c("bipolar", "mdd", "schizo", "opioids", "suicide")))
  
  df_y_res %>% 
    ggplot(aes(x = covariate, y = weight)) +
    geom_point(size = 3, shape = 1) +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    ylim(c(-0.075, 0.075)) +
    xlab("") +
    ggtitle("Y coefficient weights")
  
  # LEVEL 1 X  
  df_x_res <- tibble(ensembl_gene_id = colnames(X.mat),
                     weight = model_mat$wX) %>% 
    left_join(df_modules_filt) %>% 
    mutate(abs_weight = abs(weight)) %>% 
    arrange(-abs_weight) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = unique(.$ensembl_gene_id)))
  
  df_x_res %>% 
    #filter(abs_weight > 0.0005) %>%
    ggplot(aes(x = ensembl_gene_id, y = abs_weight, color = module)) +
    geom_point() +
    scale_color_manual(values = c(df_x_res %>% 
                                    #head(100) %>% 
                                    dplyr::count(module, color) %>% 
                                    dplyr::select(-n) %>% 
                                    pull(color))
    ) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("absolute value of X coefficient weights")
  
  df_x_res %>% 
    slice_max(abs_weight, n = 20) %>% 
    arrange(weight) %>% mutate(gene_symbol = factor(gene_symbol, levels = .$gene_symbol)) %>% 
    ggplot(aes(x = weight, y = gene_symbol, fill = weight)) +
    geom_col() +
    #coord_flip() +
    scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC") +
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
  
  df_x_res %>% 
    filter(abs_weight > 0.0005) %>% 
    count(module) %>% 
    left_join(df_modules_filt %>% 
                count(module) %>% 
                dplyr::rename("module_n" = "n")
    ) %>% 
    mutate(perc_mod = n/module_n * 100) %>% 
    arrange(-perc_mod)
  
  
  
  
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
  # X.mat <- X.mat[x_group != "M0"]
  # x_group <- x_group[x_group != "M0"]
  
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
  
  save(l_grcca_train, file = paste0(base_dir, "objects/", prefix2, "_grcca_train.RDS")) 
  
  df_train_res <- map(l_grcca_train, 4) %>% 
    map(2) %>% 
    bind_rows %>% 
    bind_cols(df_hyper, .) %>% 
    pivot_longer(contains("comp"), names_to = "comp", values_to = "r")
  
  
# RUN ACROSS VALIDATION DATA  
  # load("objects/20230206_grcca_train.RDS")
  
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
  
  save(l_grcca_test, file = paste0(base_dir, "objects/", prefix2, "_grcca_test.RDS"))

  
  df_test_res <- map(l_grcca_test, 4) %>% 
    bind_rows  
  
  

# plot parameters & variate correlations ----------------------------------

  
# LOAD DATA
  
  # train res
  # load("objects/20230215_grcca_train.RDS") # l_grcca_train

  df_train_res <- map(l_grcca_train, 4) %>%
    map(2) %>%
    bind_rows %>%
    bind_cols(df_hyper %>% filter(fold %in% 1:7), .) %>%
    pivot_longer(contains("comp"), names_to = "comp", values_to = "r")

  rm(list = "l_grcca_train")

  # test res
  # load("objects/20230215_grcca_test.RDS") # l_grcca_test

  df_test_res <- map(l_grcca_test, 4) %>%
    bind_rows

  rm(list = "l_grcca_test")



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
                 lambda1 = 100000, 
                 lambda2 = 0,
                 mu1 = 0.01,
                 mu2 = 0)
  
# CHECK THAT IT WORKED
  as.matrix(x_training) %*% grcca$x.coefs %>% head
  grcca$x.vars %>% head
  
  as.matrix(y_training) %*% grcca$y.coefs %>% head
  grcca$y.vars %>% head
  
  cor(as.matrix(x_training) %*% grcca$x.coefs, as.matrix(y_training) %*% grcca$y.coefs)
  grcca$cors
  
# CALCULATE COEFFICIENTS ON TEST  
  cor(as.matrix(x_training) %*% grcca$x.coefs, as.matrix(y_training) %*% grcca$y.coefs)
  # 0.57!!
  
# explore results ---------------------------------------------------------

    
  
# EXPLORE COEFFICIENT VECTORS  
  df_y_coefs <- grcca$y.coefs %>% 
    as.data.frame %>% 
    rownames_to_column("dx") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    pivot_longer(starts_with("can"), names_to = "can_comp", values_to = "coef") %>% 
    mutate(can_comp = str_remove(can_comp, "can_comp")) %>% 
    mutate(dx = factor(dx, levels = c("bipolar", "mdd", "schizo", "age_death", "anti_anxiety")))
  
  df_x_coefs <- grcca$x.coefs %>% 
    as.data.frame %>% 
    rownames_to_column("ensembl_transcript_id") %>% 
    as_tibble() %>% 
    pivot_longer(starts_with("can"), names_to = "can_comp", values_to = "coef") %>% 
    mutate(can_comp = str_remove(can_comp, "can.comp")) %>% 
    
    left_join(df_grcca_x %>% 
                dplyr::select(ensembl_transcript_id, module),
              by = "ensembl_transcript_id"
    ) %>% 
    distinct() 
  
  df_y_coefs %>% 
    filter(can_comp == 1) %>% 
    
    ggplot(aes(x = dx, y = coef, color = dx)) +
    geom_point(size = 3) +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    ggtitle("y coefficients of canonical component 1")
    facet_wrap( ~ can_comp, scales = "free")
  
# PLOT X COEFFICIENTS TOGETHER
  
    # assign color vector for modules
    col_vector <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"), "midnightblue", "red", "cyan", "magenta")
    col_vector <- col_vector[!col_vector %in% c("#B3B3B3", "#666666")]
    
    df_x_coefs %>% 
      
      filter(module != 0) %>% 
      
      filter(can_comp == 1) %>% 
      arrange(module, coef) %>% 
      mutate(ensembl_transcript_id = factor(ensembl_transcript_id, levels = .$ensembl_transcript_id)) %>% 
      
      ggplot(aes(x = ensembl_transcript_id, y = coef, color = module)) +
      geom_point(shape = 1) +
      geom_hline(yintercept = 0, lty = 2, color = "black") +
      scale_fill_manual(values = col_vector) +
      scale_color_manual(values = col_vector) +
      theme(axis.text.x = element_blank()) +
      ggtitle("Transcript weights for canonical component 1")
    
    
  
# TOP MODULES
    
    # plot
    df_x_coefs %>% 
      filter(can_comp == 1) %>% 
      group_by(module) %>% 
      mutate(mean_coef = mean(coef),
             sd_coef = sd(coef)) %>% 
      arrange(-abs(mean_coef), coef) %>% 
      
      left_join(df_transcript_gene_go_term) %>% 
      mutate(gene_symbol = ifelse(is.na(external_gene_name), 
                                  ensembl_transcript_id, 
                                  external_gene_name)) %>% 
      
      mutate(ensembl_transcript_id = factor(ensembl_transcript_id, levels = unique(.$ensembl_transcript_id)),
             module = factor(module, levels = unique(.$module))) %>% 
      mutate(label = ifelse(abs(coef) > mean_coef + 3*sd_coef, gene_symbol, NA)) %>% 
      dplyr::select(ensembl_transcript_id, module, coef, mean_coef, sd_coef, label) %>% 
      distinct() %>% 
      
      filter(module != 0) %>% 
      
      ggplot(aes(x = ensembl_transcript_id, y = coef)) +
      geom_point(aes(color = module), shape = 1) +
      
      # geom_text_repel(aes(label = label)) +
      
      geom_hline(aes(yintercept = mean_coef), color = "red") +
      geom_hline(aes(yintercept = mean_coef + 3*sd_coef), lty = 2, color = "black") +
      geom_hline(aes(yintercept = mean_coef - 3*sd_coef), lty = 2, color = "black") +
      
      scale_fill_manual(values = col_vector) +
      scale_color_manual(values = col_vector) +
      
      facet_wrap(vars(module), scales = "free") +
      
      theme(axis.text.x = element_blank()) +
      ggtitle("Transcript weights arrange by abs(mean_coef)")
    

    
# Z-SCORE - IDENTIFY TOP GENES & MODULES
    df_top_genes <- df_x_coefs %>% 
      filter(can_comp == 1) %>% 
      mutate(mean_coef = mean(coef),
             sd_coef = sd(coef),
             outlier = ifelse(abs(coef) > mean_coef + 3*sd_coef, 1, 0)) %>% 
      filter(outlier == 1) %>% 
      left_join(df_transcript_gene_go_term) %>% 
      mutate(gene_symbol = ifelse(is.na(external_gene_name), 
                                  ensembl_transcript_id, 
                                  external_gene_name)) %>% 
      dplyr::select(gene_symbol, module, coef) %>% 
      distinct()
    
    df_top_genes %>% 
      count(module) %>% 
      left_join(
        
        df_modules_filt %>% 
          count(module) %>% 
          dplyr::rename("module_n" = "n")
        
      ) %>% 
      mutate(perc_of_mod = n/module_n * 100) %>% 
      arrange(-perc_of_mod)
    
    # plot
    df_top_genes %>% 
      group_by(module) %>% 
      mutate(mean_coef = mean(coef)) %>% 
      arrange(-abs(mean_coef), coef) %>%
      mutate(module = factor(module, unique(.$module)),
             gene_symbol = factor(gene_symbol, unique(.$gene_symbol))) %>% 

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
      filter(module == 14) %>% 
      arrange(-abs(coef))
    
    
    
# DRUGGABLE TARGETS
    df_scz_genome <- readxl::read_excel("data/eva/lago_bahn_2022_scz_genome.xlsx") %>% 
      row_to_names(1) %>% 
      clean_names
    
    df_scz_genome$risk_gene_hgnc %>% length
    overlap_genes <- intersect(df_top_genes$gene_symbol, df_scz_genome$risk_gene_hgnc)
    
    df_top_genes %>% filter(gene_symbol %in% overlap_genes)
    
    phyper(q = length(overlap_genes), 
           m = ncol(X.mat),
           n = length(df_scz_genome$risk_gene_hgnc),
           k = length(df_top_genes$gene_symbol),
           lower.tail = TRUE)
    
    
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
  
