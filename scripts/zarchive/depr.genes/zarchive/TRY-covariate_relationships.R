

# libraries ---------------------------------------------------------------

  library(tidyverse)
  library(patchwork)
  library(janitor)
  library(RColorBrewer)
  library(tidymodels)
  library(tidytext)
  library(tidymodels)
  library(vip)
  library(workflowsets)
  library(ggpubr)
  library(corrr)
  library(pheatmap)
  library(cluster)




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
  load(paste0(base_dir, "data/covariates/185_all_covariates_clean.Rdata"))




# covariate correlations --------------------------------------------------


  technical_covariates <- c("mapped_percent", "five_prime_three_prime_bias", "rin_acsg", 
                            "rn_aextraction_batch", "library_batch", "pmi", "g_cpercent", 
                            "ph")
  
  confounds <- c("age_death", "brain_weight", "height", "weight")

  # remove manner because it's same as suicide
  # remove BMI because it's same as weight
  
# CALCULATE CORRELATIONS  
  m_covariates_cor <- df_covariates_all_clean_numeric %>% 
    dplyr::select(-all_of(technical_covariates), -manner, -bmi) %>%
    column_to_rownames("sample") %>% 
    cor(use = "pairwise.complete.obs")
  diag(m_covariates_cor) <- 0

# PLOT HEATMAP
  pheatmap(m_covariates_cor,
           color = rev(brewer.pal(9, "RdBu")))

# COVARIATES WITH SIGNIFICANT CORRELATION
  df_covariate_correlations <- m_covariates_cor %>% 
    as.data.frame %>% 
    rownames_to_column("term1") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), names_to = "term2", values_to = "r") %>% 
    
    # remove duplicates in reverse
    group_by(grp = paste0(pmax(term1, term2), sep = "_", pmin(term1, term2))) %>% 
    slice(1) %>% 
    ungroup %>% 
    dplyr::select(-grp) %>% 
    
    # get absolute value of correlation
    mutate(abs_r = abs(r)) %>% 
    arrange(-abs_r) %>% 
    
    # calculate z_score
    mutate(z_score = (r - mean(r))/sd(r))
    
# PLOT HISTOGRAM  
  df_covariate_correlations %>% 
    filter(abs(z_score) > 2) %>% 
    arrange(abs(z_score))
  
  df_covariate_correlations %>% 
    ggplot(aes(x = r)) + 
    geom_histogram() +
    
    geom_vline(aes(xintercept = 0.385), lty = 2, color = "red") +
    geom_vline(aes(xintercept = -0.33), lty = 2, color = "red") +
    annotate(geom = "text", label = "z-score > 2 →", x = 0.45, y = 75, hjust = 0) +
    annotate(geom = "text", label = "← z-score < -2", x = -0.4, y = 75, hjust = 1) +
    
    xlim(c(-1, 1)) +
    ylab("") +
    ggtitle("histogram of covariate correlations")
  
# PLOT RELATIONSHIPS BETWEEN SIGNIFICANT COVARIATES 
  df_covariate_correlations %>% 
    filter(abs(z_score) > 2) %>% 
    arrange(abs(z_score))

  

# covariate PCA ---------------------------------------------------------------------


# CLUSTER COVARIATES
  
  # silhouette
  dist <- dist(m_covariates_cor)
  sil_width <- map_dbl(2:20, function(k) {
    model <- pam(dist, k = k)
    model$silinfo$avg.width
  })
  df_sil <- tibble(k = 2:20, sil_width = sil_width)
  df_sil %>% ggplot(aes(x = k, y = sil_width)) +
    geom_point() +
    geom_line() +
    #ylim(c(0.3, 0.6)) +
    ggtitle("Silhouette widths")
  
  df_clust_covar <- cutree(hclust(dist), k = 9) %>% 
    enframe() %>% 
    dplyr::rename("term" = "name", "cluster" = "value") %>% 
    mutate(cluster = factor(cluster))
  
# RUN PCA  
  pca <- df_covariates_all_clean_numeric %>% 
    dplyr::select(-all_of(technical_covariates), -manner, -bmi) %>%
    mutate_if(is.numeric, ~ ifelse(is.na(.x), 0, .x)) %>% 
    column_to_rownames("sample") %>% 
    prcomp(center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("term") %>% 
    as_tibble() %>% 
    
    left_join(df_clust_covar)
    
  df_pca %>% 
    ggplot(aes(x = PC1, y = PC2, color = cluster)) +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_text_repel(aes(label = term)) +
    ggtitle("covariate PCA")
    
  
# Drug correlations --------------------------------------------------
  
  
  technical_covariates <- c("mapped_percent", "five_prime_three_prime_bias", "rin_acsg", 
                            "rn_aextraction_batch", "library_batch", "pmi", "g_cpercent", 
                            "ph")
  
  confounds <- c("age_death", "brain_weight", "height", "weight")
  
  drugs <- c("nicotine", "cocaine", "alcohol", "anti_anxiety", "opioids",
             "cannabis", "major_stimulants", "anti_depressant", "antipsychotics",
             "mood_stabilizers", "anti_epileptics", "anti_histamines", "anticholinergics",
             "smoker", "other_psychotropic_drug", "non_psychiatric")
  
  # CALCULATE CORRELATIONS  
  m_drugs_cor <- df_covariates_all_clean_numeric %>% 
    dplyr::select(sample, all_of(drugs)) %>%
    column_to_rownames("sample") %>% 
    cor(use = "pairwise.complete.obs")
  diag(m_drugs_cor) <- 0
  
  # PLOT HEATMAP
  pheatmap(m_drugs_cor,
           color = rev(brewer.pal(9, "RdBu")))
  
  # COVARIATES WITH SIGNIFICANT CORRELATION
  df_drug_correlations <- m_drugs_cor %>% 
    as.data.frame %>% 
    rownames_to_column("term1") %>% 
    as_tibble() %>% 
    pivot_longer(2:ncol(.), names_to = "term2", values_to = "r") %>% 
    
    # remove duplicates in reverse
    group_by(grp = paste0(pmax(term1, term2), sep = "_", pmin(term1, term2))) %>% 
    slice(1) %>% 
    ungroup %>% 
    dplyr::select(-grp) %>% 
    
    # get absolute value of correlation
    mutate(abs_r = abs(r)) %>% 
    arrange(-abs_r) %>% 
    
    # calculate z_score
    mutate(z_score = (r - mean(r))/sd(r))
  
  # PLOT HISTOGRAM  
  df_drug_correlations %>% 
    filter(abs(z_score) > 2) %>% 
    arrange(abs(z_score))
  
  df_drug_correlations %>% 
    ggplot(aes(x = r)) + 
    geom_histogram() +
    
    geom_vline(aes(xintercept = 0.55), lty = 2, color = "red") +
    annotate(geom = "text", label = "z-score > 2 →", x = 0.6, y = 25, hjust = 0) +

    xlim(c(-1, 1)) +
    ylab("") +
    ggtitle("histogram of drug correlations")
  
# PLOT RELATIONSHIPS BETWEEN SIGNIFICANT COVARIATES 
  df_drug_correlations %>% 
    filter(abs(z_score) > 2) %>% 
    arrange(abs(z_score))
  
  # drug class 1: anti_histamines, non_psychiatric, other_psychotropic_drug
  # drug class 2: anti_anxiety, mood_stabilizers, anti_epileptics
  # drug class 3: smoker, nicotine

  
# covariate PCA ---------------------------------------------------------------------
  
  
# CLUSTER DRUGS
  
  # silhouette
  dist <- dist(m_drugs_cor)
  sil_width <- map_dbl(2:10, function(k) {
    model <- pam(dist, k = k)
    model$silinfo$avg.width
  })
  df_sil <- tibble(k = 2:10, sil_width = sil_width)
  df_sil %>% ggplot(aes(x = k, y = sil_width)) +
    geom_point() +
    geom_line() +
    ggtitle("Silhouette widths")
  
  df_clust_drugs <- cutree(hclust(dist), k = 6) %>% 
    enframe() %>% 
    dplyr::rename("term" = "name", "cluster" = "value") %>% 
    mutate(cluster = factor(cluster))
  
# RUN PCA  
  pca <- df_covariates_all_clean_numeric %>% 
    dplyr::select(sample, all_of(drugs)) %>%
    mutate_if(is.numeric, ~ ifelse(is.na(.x), 0, .x)) %>% 
    column_to_rownames("sample") %>% 
    prcomp(center = TRUE, scale = TRUE)
  
  df_pca <- pca$rotation %>% 
    as.data.frame %>% 
    rownames_to_column("term") %>% 
    as_tibble() %>% 
    
    left_join(df_clust_drugs) %>% 
    
    mutate(drug_class = case_when(
      
      term %in% c("anti_histamines", "non_psychiatric", "other_psychotropic_drug") ~ "other",
      term %in% c("anti_anxiety", "mood_stabilizers", "anti_epileptics") ~ "mood_meds",
      term %in% c("smoker", "nicotine", "alcohol", "cannabis", "cocaine", "opioids", "major_stimulants") ~ "drugs_of_abuse",
      term %in% c("anticholinergics", "anti_depressant", "antipsychotics" ) ~ "other_psych_meds"
      
    ))
  
  df_pca %>% 
    ggplot(aes(x = PC1, y = PC2, color = drug_class)) +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_text_repel(aes(label = term)) +
    ggtitle("covariate PCA")
  
  
  
  
  
  
  
  