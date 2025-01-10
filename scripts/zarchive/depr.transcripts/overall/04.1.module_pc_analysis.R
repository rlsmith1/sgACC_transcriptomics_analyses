
########################################################################################

# OVERALL TRANSCRIPTS: Associate module PCs with covariates 

########################################################################################


# libraries ---------------------------------------------------------------

  library(tidyverse)
  library(patchwork)
  library(janitor)
  library(RColorBrewer)
  library(tidymodels)
  library(tidytext)
  library(R.matlab)


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
prefix <- "20230301_185samples_70ktranscripts_vst_qSVA123567_MP_RNAbatch_Race_resids"
soft_power <- 3

# LOAD OBJECTS  
  load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
  load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

# IDENTIFY MODULE SET OF INTEREST  
  df_modules_filt <- df_modules %>% 
    filter(min_size == 40 & cut_height == 0.996) %>% 
    dplyr::select(ensembl_transcript_id, module)

# COVARIATES
  load("data/covariates/185_all_covariates_clean.Rdata")
  


# run PCA on each module --------------------------------------------------


# RUN PCA FOR EACH MODULE
  df_mods_pca <- df_modules_filt %>% 
    left_join(df_vsd_regress) %>% 
    pivot_wider(id_cols = c(ensembl_transcript_id, module), 
                names_from = sample, 
                values_from = resids) %>% 
    group_by(module) %>% 
    nest() %>% 
    
    mutate(pca = map(.x = data,
                     .f = ~ .x %>% 
                       as.data.frame %>%
                       column_to_rownames("ensembl_transcript_id") %>% 
                       prcomp(center = TRUE, scale = TRUE)
    )) %>% 
    
    mutate(sample_pcs = map(.x = pca,
                            .f = ~.x$rotation %>% 
                              as.data.frame %>% 
                              rownames_to_column("sample") %>% 
                              as_tibble() %>% 
                              clean_names() %>% 
                              pivot_longer(matches("*[0-9]"), 
                                           names_to = "pc", 
                                           values_to = "pc_score") %>% 
                              
                              left_join(
                                
                                summary(.x) %>% 
                                  .$importance %>% 
                                  as.data.frame %>% 
                                  rownames_to_column("metric") %>% 
                                  as_tibble() %>% 
                                  pivot_longer(2:ncol(.), 
                                               values_to = "value", 
                                               names_to = "pc") %>% 
                                  mutate(pc = tolower(pc)) %>% 
                                  pivot_wider(id_cols = pc, 
                                              names_from = metric, 
                                              values_from = value) %>% 
                                  clean_names()
                                
                              ) %>% 
                              filter(proportion_of_variance > 0.01)
    )
    ) %>% 
    arrange(module) %>% 
    unnest(cols = c(sample_pcs)) %>% 
    dplyr::select(module, sample, pc, proportion_of_variance, pc_score) %>% 
    
    left_join(df_covariates_all_clean_numeric %>% 
                pivot_longer(2:ncol(.), 
                             names_to = "covariate", 
                             values_to = "covariate_val"),
              by = "sample")
  
  
  
# quantify module-trait associations --------------------------------------
  
  
# CORRELATE WITH COVARIATES
  df_mods_pc_covar_cor <- df_mods_pca %>% 
    
    group_by(module, pc, covariate) %>% 
    nest() %>% 
    
    mutate(pearsons_r = map(.x = data,
                            .f = ~ cor.test(.x$pc_score, .x$covariate_val)$estimate
    ),
    p_val = map(.x = data,
                .f = ~ cor.test(.x$pc_score, .x$covariate_val)$p.value
    )) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    ungroup %>% 
    mutate(q_val = p.adjust(p_val, method = "BH")) %>% 
    mutate(pc = str_remove(pc, "pc") %>% 
             as.numeric %>% 
             as.factor)
  
  df_mods_pc_covar_cor %>% 
    filter(q_val < 0.05) %>% 
    arrange(pc, module) %>% 
    print(n=nrow(.))
  
# PLOT
  df_mods_pc_covar_cor %>% 
    
    unnest(cols = c(data)) %>% 
    group_by(module) %>% 
    filter(proportion_of_variance > 0.02) %>% # filter(q_val < 0.05) %>% select(module, pc, covariate) %>% distinct

    mutate(pearsons_r = ifelse(p_val > 0.05, NA, pearsons_r)) %>% 
    mutate(labels = case_when(
      
      q_val < 0.05 ~ "*",
      q_val < 0.01 ~ "**",
      q_val < 0.001 ~ "***"
      
      
    )) %>% 
    
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = pearsons_r)) +
    geom_text(aes(label = labels), color = "black", size = 4) +
    scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC", na.value = "white",
                         limits = c(-1, 1)) +
    facet_wrap(~ module, scales = "free") +
    labs(x = "PC explained var > 2%") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 8),
          strip.text = element_text(size = 10)) +
    ggtitle("Module PC-covariate correlations")
  


  
# module eigengene-diagnosis correlation ----------------------------------


# GET EIGENGENES
  df_kme <- df_mods_pca %>% 
    filter(pc == "pc1") %>% 
    pivot_wider(id_cols = c(module, sample, pc_score), 
                names_from = covariate, 
                values_from = covariate_val) %>% 
    dplyr::rename("kme" = "pc_score") %>% 
    mutate(dx = case_when(
      
      multi_nomial_dx == 1 ~ "control",
      multi_nomial_dx == 2 ~ "bipolar",
      multi_nomial_dx == 3 ~ "mdd",
      multi_nomial_dx == 4 ~ "schizo"
      
    ) %>% 
      factor(levels = c("control", "bipolar", "mdd", "schizo")),
    .before = 4)

# PLOT EIGENGENES ACROSS DX
  df_kme %>% 
    
    ggplot() +
    geom_boxplot(aes(x = kme, y = -2.5, color = dx), width = 5) +
    geom_density(aes(x = kme, fill = dx), alpha = 0.5) +
    facet_wrap(~module) +
    labs(y = "") +
    ggtitle("module eigengene distributions by diagnosis in each module")
  
  

# LINEAR REGRESSION OF MES AND DIAGNOSIS

  df_lm_dx <- df_kme %>% 
    group_by(module) %>% 
    nest() %>% 
    mutate(
      
      lm = map(
        
        .x = data,
        .f = ~ lm(kme ~ dx, data = .x)
        
      ),
      
      lm_sum = map(
        
        .x = lm,
        .f = ~ summary(.x)
        
      ),
      
      lm_res = map(
        
        .x = lm,
        .f = ~ tidy(.x)
        
      )
      
    ) %>% 
    unnest(cols = c(lm_res)) %>% 
    clean_names %>% 
    filter(term != "(Intercept)") %>% 
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           term = str_remove(term, "dx"))
  
  df_lm_dx %>% 
    filter(p_value < 0.05) %>% 
    dplyr::select(module, term, estimate, p_value, p_adj)



    
# select best model per module --------------------------------------------
  
  library(tidymodels)
  library(vip)
  library(workflowsets)
  
  
# TIDYMODELS COMPARISON
  
  # recipes to test on all modules
  rec1 <- recipe(kme ~ dx, data = df_kme)
  rec2 <- recipe(kme ~ dx + opioids, data = df_kme)
  rec3 <- recipe(kme ~ dx + nicotine, data = df_kme)
  rec4 <- recipe(kme ~ dx + nicotine + opioids, data = df_kme)
  rec5 <- recipe(kme ~ dx + brain_weight, data = df_kme)
  rec6 <- recipe(kme ~ dx + brain_weight + opioids, data = df_kme)
  rec7 <- recipe(kme ~ dx + age_death, data = df_kme)
  rec8 <- recipe(kme ~ dx + brain_weight + age_death, data = df_kme)
  rec9 <- recipe(kme ~ dx + brain_weight + suicide, data = df_kme)
  rec10 <- recipe(kme ~ dx + opioids + suicide, data = df_kme)
  rec11 <- recipe(kme ~ dx + suicide, data = df_kme)
  rec12 <- recipe(kme ~ dx + suicide + opioids + brain_weight, data = df_kme)
  rec13 <- recipe(kme ~ dx + age_death + opioids, data = df_kme)
  
  # run linear models only
  lm_model <- linear_reg() %>% 
    set_engine('lm') %>% # adds lm implementation of linear regression
    set_mode('regression')
  
  # workflow set
  all_models <-  workflow_set(
    preproc = list(rec1 = rec1, 
                   rec2 = rec2, 
                   rec3 = rec3,
                   rec4 = rec4,
                   rec5 = rec5,
                   rec6 = rec6,
                   rec7 = rec7,
                   rec8 = rec8,
                   rec9 = rec9,
                   rec10 = rec10,
                   rec11 = rec11,
                   rec12 = rec12,
                   rec13 = rec13
    ),
    models = list(lm = lm_model),
    cross = TRUE
  )
  
  
  df_mods_best_rec <- map_dfr(
    
    .x = unique(df_kme$module),
    .f = ~ list( {
      
      # select module
      df_mod <- df_kme %>% 
        filter(module == .x)
      
      # split data
      set.seed(123)
      mod_split <- initial_split(df_mod, prop = 0.75, strata = dx)
      mod_training <- training(mod_split)
      mod_testing <- testing(mod_split)
      
      # fit each workflow to data
      all_models_fit <- all_models %>% 
        mutate(fit = 
                 
                 map(
                   
                   .x = info,
                   .f = ~ .x$workflow %>% 
                     .[[1]] %>% 
                     last_fit(split = mod_split)
                   
                 )
               
        )
      
      # Obtain performance metrics and predictions on test data
      df_mod_models_res <- all_models_fit %>% 
        mutate(metrics = map(
          
          .x = fit,
          .f = ~ .x %>% collect_metrics
          
        ),
        
        predictions = map(
          
          .x = fit,
          .f = ~ .x %>% collect_predictions
          
        )
        
        ) %>% 
        
        unnest(cols = c(metrics)) %>% 
        dplyr::select(-.config) %>% 
        unnest(cols = c(predictions))
      
      best_rec <- df_mod_models_res %>% 
        dplyr::select(wflow_id, .metric, .estimate) %>% 
        distinct %>% 
        filter(.metric == "rsq") %>% 
        arrange(-.estimate) %>% 
        head(1) %>% pull(wflow_id) %>% str_remove("_lm")
      
      return(
        
        tibble(
          
          mod = .x,
          best_rec = best_rec
          
        )
        
      )
      
    })
    
  )  
  
  df_mods_best_rec %>% 
    dplyr::count(best_rec) %>% 
    arrange(-n)
  
  # final fit with best recipe
  df_res <- tibble()
  for (i in 1:nrow(df_mods_best_rec)) {
    
    df_current_best_rec <- df_mods_best_rec[i,]
    
    if (df_current_best_rec$best_rec == "rec1") {
      
      best_rec <- rec1
      
    } else if (df_current_best_rec$best_rec == "rec2") {
      
      best_rec <- rec2
      
    } else if (df_current_best_rec$best_rec == "rec3") {
      
      best_rec <- rec3
      
    } else if (df_current_best_rec$best_rec == "rec4") {
      
      best_rec <- rec4
      
    } else if (df_current_best_rec$best_rec == "rec5") {
      
      best_rec <- rec5
      
    } else if (df_current_best_rec$best_rec == "rec6") {
      
      best_rec <- rec6
      
    } else if (df_current_best_rec$best_rec == "rec7") {
      
      best_rec <- rec7
      
    }  else if (df_current_best_rec$best_rec == "rec8") {
      
      best_rec <- rec8
      
    }  else if (df_current_best_rec$best_rec == "rec9") {
      
      best_rec <- rec9
      
    }  else if (df_current_best_rec$best_rec == "rec10") {
      
      best_rec <- rec10
      
    }  else if (df_current_best_rec$best_rec == "rec11") {
      
      best_rec <- rec11
      
    }  else if (df_current_best_rec$best_rec == "rec12") {
      
      best_rec <- rec12
      
    }  else if (df_current_best_rec$best_rec == "rec13") {
      
      best_rec <- rec13
      
    }
    
    mod_wflow <- workflow() %>% 
      add_recipe(best_rec) %>% 
      add_model(lm_model)
    
    df_tmp <- mod_wflow %>% 
      fit(data = df_kme %>% filter(module == df_current_best_rec$mod)) %>% 
      tidy %>% 
      filter(term != "(Intercept)") %>% 
      mutate(module = df_current_best_rec$mod, 
             best_rec = df_current_best_rec$best_rec,
             .before = 1) %>% 
      clean_names()
    
    df_res <- df_res %>% 
      bind_rows(df_tmp)
    
  }
  
  df_res %>% 
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>% 
    distinct() %>% 
    filter(p_value < 0.05) #%>% 
  write.csv("outputs/20230308_univariate_lm_bestRec.csv",
            row.names = FALSE)
  
  
  # opioids?
  df_kme %>% 
    filter(module %in% c(1, 8, 10, 15, 16) & !is.na(opioids)) %>% 
    mutate(opioids = ifelse(opioids == 1, "yes", "no")) %>% 
    ggplot(aes(x = dx, y = kme, color = opioids)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                               dodge.width = 0.75)) +
    facet_wrap(vars(module), nrow = 3) +
    labs(x = "diagnosis") +
    ggtitle(aes("Module eigengene by diagnosis and opioid use"))
  
  
  
  
  
  
  
  
  


  

