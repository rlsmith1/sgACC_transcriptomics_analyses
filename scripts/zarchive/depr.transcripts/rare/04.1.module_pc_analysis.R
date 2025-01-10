
# !!!!!! RARE !!!!!!

########################################################################################

# Associate module PCs with covariates 

########################################################################################


# libraries ---------------------------------------------------------------

  library(tidyverse)
  library(patchwork)
  library(janitor)
  library(RColorBrewer)
  library(tidymodels)
  library(tidytext)
  library(R.matlab)
  library(ggpubr)


# set theme for plots -----------------------------------------------------

  theme_set(theme_bw() +
              theme(plot.title = element_text(size = 18),
                    axis.title = element_text(size = 15),
                    axis.text = element_text(size = 12),
                    strip.text = element_text(size = 15),
                    legend.title = element_text(size = 15),
                    legend.text = element_text(size = 12)))
  



# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "15Nov2023_RARE_regress_qSV1.7_ageDeath_RNAbatch_GCperc"

soft_power <- 2
minimum_size <- 35
tree_cut_height <- 0.9936

# LOAD OBJECTS 
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

# IDENTIFY MODULE SET OF INTEREST  
df_modules_filt <- df_modules %>% 
  filter(min_size == minimum_size & cut_height == tree_cut_height) %>% 
  unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
  arrange(mod_set, module)

# COVARIATES
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_clean

# DRUG MCA RESULTS
load(paste0(base_dir, "objects/drug_MCA_results.Rdata"))

# TRANSCRIPT SYMBOL TO ENSEMBL TRANSCRIPT ID MAPPING
load(paste0(base_dir, "objects/df_hsapiens_genome.RDS")) # df_hsapiens_genome

df_ensembl_to_symbol <- df_hsapiens_genome %>% 
  dplyr::select(transcript_id, transcript_name) %>% 
  distinct() %>% 
  dplyr::rename("transcript_symbol" = "transcript_name", "ensembl_transcript_id" = "transcript_id")

df_gene_to_transcript <- df_hsapiens_genome %>% 
  dplyr::select(transcript_id, transcript_name, gene_id, gene_name) %>% 
  distinct() %>% 
  dplyr::rename("gene_symbol" = "gene_name", "transcript_symbol" = "transcript_name", 
                "ensembl_gene_id" = "gene_id", "ensembl_transcript_id" = "transcript_id"
  ) %>% 
  filter(!is.na(ensembl_transcript_id))

df_modules_filt <- df_modules_filt %>% 
  left_join(df_gene_to_transcript) %>% 
  dplyr::select(ensembl_transcript_id, ensembl_gene_id, gene_symbol, module, color) %>% 
  distinct()




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
    
    ggplot(aes(x = kme, y = dx, color = dx)) +
    geom_violin() +
    stat_summary(fun = "median", geom = "crossbar", width = 1) +
    geom_point(position = position_jitter(height = 0.25), shape = 1, alpha = 0.7) +
    facet_wrap(vars(module)) +
    labs(y = "") +
    ggtitle("module eigengene distributions by diagnosis in each module") +
    theme(strip.text = element_text(size = 10),
          legend.position = "none")  
    

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
    filter(p_adj < 0.05) %>% 
    dplyr::select(module, term, estimate, p_value, p_adj) %>% 
    
    write.csv("outputs/20230316_RareTx_univariate_lm_noCovar.csv",
              row.names = FALSE)
  

  
  
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
    unnest(cols = c(data)) %>% 
    filter(proportion_of_variance > 0.02) %>% 
    dplyr::select(module, pc, proportion_of_variance, covariate, pearsons_r, p_val) %>% 
    distinct() %>% 
    filter(p_val < 0.05) %>% 
    arrange(pc, module) %>% 
    count(covariate) %>% 
    filter(n > 1) %>% 
    arrange(-n) %>% 
    filter(!(covariate %in% c("g_cpercent", "ph", "pmi", "rin_acsg")))

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
  

  

# select best model per module using covariates --------------------------------------------
  
  
# IDENTIFY POTENTIAL COVARIATES BASED ON PC CORRELATION
  covariates <- df_mods_pc_covar_cor %>% 
    unnest(cols = c(data)) %>% 
    filter(proportion_of_variance > 0.02) %>% 
    dplyr::select(module, pc, proportion_of_variance, covariate, pearsons_r, p_val) %>% 
    distinct() %>% 
    filter(p_val < 0.05) %>% 
    arrange(pc, module) %>% 
    count(covariate) %>% 
    filter(n > 1) %>% 
    arrange(-n) %>% 
    filter(covariate != "multi_nomial_dx") %>% 
    
    # remove technical covariates
    filter(!(covariate %in% c("g_cpercent", "ph", "pmi", "rin_acsg"))) %>% 
    
    # remove antipsychotics bc confounded with SCZ,
    # manner of death because suicide accounts for that,
    # and bmi because it correlates with weight
    filter(!(covariate %in% c("antipsychotics", "manner", "bmi"))) %>% 
    pull(covariate)
  
# ROUND 1: DX PLUS ONE COVARIATE
  
  recs <- map(.x = covariates, 
              .f = ~ recipe(as.formula(paste0("kme ~ dx + ", .x)), data = df_kme))
  
  # run linear models only
  lm_model <- linear_reg() %>% 
    set_engine('lm') %>% # adds lm implementation of linear regression
    set_mode('regression')
  
  # workflow set
  all_models <-  workflow_set(
    preproc = recs,
    models = list(lm = lm_model),
    cross = TRUE
  )
  
  # select best recipe for each module based on adj r-squared on test data
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
  
  # final fit with best recipe for each module
  df_round1_res <- df_mods_best_rec %>% 
    mutate(
      
      best_rec = map(
        
        .x = best_rec,
        .f = ~ .x %>% 
          str_remove("recipe_") %>% 
          as.numeric %>% 
          recs[[.]]
        
      ),
      
      wflow = map(
        
        .x = best_rec,
        .f = ~ workflow() %>% 
          add_recipe(.x) %>% 
          add_model(lm_model)
        
      ),
      
      res = map2(
        
        .x = wflow,
        .y = mod,
        .f = ~ .x %>% 
          fit(data = df_kme %>% filter(module == .y)) %>% 
          tidy %>% 
          clean_names %>% 
          filter(term != "(Intercept)")
        
      )
      
    ) %>% 
    unnest(cols = c(res)) %>% 
    mutate(p_adj = p.adjust(p_value, method = "BH"))
  
  # identify significant covariates
  df_round1_res %>% 
    filter(p_value < 0.05)

# ROUND 2: ALL COMBOS OF 2 SIGNIFICANT COVARIATES  
  
 covariates2 <- expand_grid(cov1 = covariates,
             cov2 = covariates) %>% 
   filter(cov1 != cov2) %>% 
   group_by(grp = paste(pmax(cov1, cov2), pmin(cov1, cov2), sep = "_")) %>%
   slice(1) %>%
   ungroup() %>% 
   unite(formula, c(cov1, cov2), sep = " + ") %>% 
   pull(formula)
  
 recs2 <- map(.x = c(covariates, covariates2), 
             .f = ~ recipe(as.formula(paste0("kme ~ dx + ", .x)), data = df_kme))
 
 # run linear models only
 lm_model <- linear_reg() %>% 
   set_engine('lm') %>% # adds lm implementation of linear regression
   set_mode('regression')
 
 # workflow set
 all_models <-  workflow_set(
   preproc = recs2,
   models = list(lm = lm_model),
   cross = TRUE
 )
 
 # select best recipe for each module based on adj r-squared on test data
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
 
 # final fit with best recipe for each module
 df_round2_res <- df_mods_best_rec %>% 
   mutate(
     
     best_rec = map(
       
       .x = best_rec,
       .f = ~ .x %>% 
         str_remove("recipe_") %>% 
         as.numeric %>% 
         recs2[[.]]
       
     ),
     
     wflow = map(
       
       .x = best_rec,
       .f = ~ workflow() %>% 
         add_recipe(.x) %>% 
         add_model(lm_model)
       
     ),
     
     res = map2(
       
       .x = wflow,
       .y = mod,
       .f = ~ .x %>% 
         fit(data = df_kme %>% filter(module == .y)) %>% 
         tidy %>% 
         clean_names %>% 
         filter(term != "(Intercept)")
       
     )
     
   ) %>% 
   unnest(cols = c(res)) %>% 
   mutate(p_adj = p.adjust(p_value, method = "BH"))
 
 # identify significant covariates
 df_round2_res %>% 
   dplyr::select(-best_rec, -wflow) %>% 
   left_join(df_mods_best_rec) %>% 
   dplyr::select(mod, best_rec, term, estimate, p_value, p_adj) %>% 
   filter(p_value < 0.05) 
  
 
# ROUND 3: ALL COMBOS OF SIGNIFICANT COVARIATES
   sig_covariates <- df_round2_res %>% 
     filter(p_value < 0.05 & !grepl("dx", term)) %>% 
     pull(term) %>% unique
   
   sig_covariates2 <- expand_grid(cov1 = sig_covariates,
                                  cov2 = sig_covariates) %>% 
     filter(cov1 != cov2) %>% 
     group_by(grp = paste(pmax(cov1, cov2), pmin(cov1, cov2), sep = "_")) %>%
     slice(1) %>%
     ungroup() %>% 
     unite(formula, c(cov1, cov2), sep = " + ") %>% 
     pull(formula)
   
   sig_covariates3 <- expand_grid(cov1 = sig_covariates,
               cov2 = sig_covariates,
               cov3 = sig_covariates) %>% 
     filter(cov1 != cov2 & cov2 != cov3 & cov1 != cov3) %>% 
     
     # remove repeats
     group_by(grp = paste(pmax(cov1, cov2, cov3), pmin(cov1, cov2, cov3), sep = "_")) %>%
     slice(1) %>%
     ungroup() %>% 

     unite(formula, c(cov1, cov2, cov3), sep = " + ") %>% 
     pull(formula)
   
   sig_covariates4 <- paste(sig_covariates, collapse = " + ")
   
   all_recs <-  map(.x = c(sig_covariates, sig_covariates2, sig_covariates3, sig_covariates4), 
                    .f = ~ recipe(as.formula(paste0("kme ~ dx + ", .x)), data = df_kme))
   
   # run linear models 
   lm_model <- linear_reg() %>% 
     set_engine('lm') %>% # adds lm implementation of linear regression
     set_mode('regression')
   
   # workflow set
   all_models <-  workflow_set(
     preproc = all_recs,
     models = list(lm = lm_model),
     cross = TRUE
   )
   
   # select best recipe for each module based on adj r-squared on test data
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
   
   # final fit with best recipe for each module
   df_round2_res <- df_mods_best_rec %>% 
     mutate(
       
       best_rec = map(
         
         .x = best_rec,
         .f = ~ .x %>% 
           str_remove("recipe_") %>% 
           as.numeric %>% 
           all_recs[[.]]
         
       ),
       
       wflow = map(
         
         .x = best_rec,
         .f = ~ workflow() %>% 
           add_recipe(.x) %>% 
           add_model(lm_model)
         
       ),
       
       res = map2(
         
         .x = wflow,
         .y = mod,
         .f = ~ .x %>% 
           fit(data = df_kme %>% filter(module == .y)) %>% 
           tidy %>% 
           clean_names %>% 
           filter(term != "(Intercept)")
         
       )
       
     ) %>% 
     unnest(cols = c(res)) %>% 
     mutate(p_adj = p.adjust(p_value, method = "BH"))
   
   # identify significant covariates
   df_round2_res %>% 
     dplyr::select(-best_rec, -wflow) %>% 
     left_join(df_mods_best_rec) %>% 
     dplyr::select(mod, best_rec, term, estimate, p_value, p_adj) %>% 
     filter(p_value < 0.05) #%>% 
   write.csv("outputs/20230316_RareTx_univariate_lm_bestRec.csv",
             row.names = FALSE)
  
      
   
   

# MODULE PLOTS ------------------------------------------------------------

    

# 5
   p_5a <- df_kme %>% 
     filter(module == 5 & !is.na(anticholinergics)) %>% 
     mutate(anticholinergics = ifelse(anticholinergics == 1, "yes", "no")) %>% 
     ggplot(aes(x = anticholinergics, y = kme, color = anticholinergics)) +
     geom_boxplot() +
     geom_point(position = position_jitterdodge(jitter.width = 0.5, 
                                                dodge.width = 0.75)) +
     
     labs(x = "anticholinergic use?") +
     theme(legend.position = "bottom")
   
   p_5b <- df_kme %>% 
     filter(module == 5 & !is.na(anticholinergics)) %>% 
     mutate(anticholinergics = ifelse(anticholinergics == 1, "yes", "no")) %>% 
     ggplot(aes(x = anticholinergics, y = kme, color = dx)) +
     geom_boxplot() +
     geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                                dodge.width = 0.75)) +
     
     labs(x = "anticholinergic use?") +
     theme(legend.position = "bottom")
   
   p_5a + p_5b + plot_annotation(title = "Module 5")
   
# 6   
   p_6a <- df_kme %>% 
     filter(module == 6 & !is.na(anti_anxiety)) %>% 
     mutate(anti_anxiety = ifelse(anti_anxiety == 1, "yes", "no")) %>% 
     ggplot(aes(x = anti_anxiety, y = kme, color = anti_anxiety)) +
     geom_boxplot() +
     geom_point(position = position_jitterdodge(jitter.width = 0.5, 
                                                dodge.width = 0.75)) +
     
     labs(x = "anti-anxiety use?") +
     theme(legend.position = "bottom")
   
   p_6b <- df_kme %>% 
     filter(module == 6 & !is.na(anti_anxiety)) %>% 
     mutate(anti_anxiety = ifelse(anti_anxiety == 1, "yes", "no")) %>% 
     ggplot(aes(x = anti_anxiety, y = kme, color = dx)) +
     geom_boxplot() +
     geom_point(position = position_jitterdodge(jitter.width = 0.5, 
                                                dodge.width = 0.75)) +
     
     labs(x = "anti-anxiety use?") +
     theme(legend.position = "bottom")
   
   p_6c <- df_kme %>% 
     filter(module == 6) %>% 
     mutate(anti_anxiety = ifelse(anti_anxiety == 1, "yes", "no")) %>% 
     
     ggplot(aes(x = age_death, y = kme)) +
     geom_point() +
     geom_smooth(method = "lm") +
     stat_cor(label.y = 0.2)
   
   p_6d <- df_kme %>% 
     filter(module == 6) %>% 
     filter(!is.na(anti_anxiety)) %>% 
     mutate(anti_anxiety = ifelse(anti_anxiety == 1, "anti-anxiety use", "no anti-anxiety use")) %>% 
     
     ggplot(aes(x = age_death, y = kme, color = dx)) +
     geom_point() +
     geom_smooth(method = "lm") +
     facet_grid(anti_anxiety ~ dx) +
     stat_cor(label.y = 0.2, size = 3) +
     theme(legend.position = "none",
           strip.text = element_text(size = 10))
   
   (p_6a + p_6c) / (p_6b + p_6d) + plot_annotation("Module 6")
  
   
# 8
   p_8a <- df_kme %>% 
     filter(module == 8 & !is.na(opioids)) %>% 
     mutate(opioids = ifelse(opioids == 1, "yes", "no")) %>% 
     ggplot(aes(x = opioids, y = kme, color = opioids)) +
     geom_boxplot(outlier.shape = NA) +
     geom_point(position = position_jitterdodge(jitter.width = 0.5, 
                                                dodge.width = 0.75)) +
     
     labs(x = "opioid use?") +
     theme(legend.position = "bottom")
   
   p_8b <- df_kme %>% 
     filter(module == 8 & !is.na(opioids)) %>% 
     mutate(opioids = ifelse(opioids == 1, "yes", "no")) %>% 
     ggplot(aes(x = opioids, y = kme, color = dx)) +
     geom_boxplot(outlier.shape = NA) +
     geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                                dodge.width = 0.75)) +
     
     labs(x = "opioid use?") +
     theme(legend.position = "bottom")
   
   p_8a + p_8b + plot_annotation(title = "Module 8")
   
   
# 12, 13, 14   
   p12_13_14 <- df_kme %>% 
     filter(module %in% c(12, 13, 14) & !is.na(anti_anxiety)) %>% 
     mutate(opioids = ifelse(opioids == 1, "yes", "no")) %>% 
     ggplot(aes(x = dx, y = kme, color = dx)) +
     geom_boxplot(outlier.shape = NA) +
     geom_point(position = position_jitterdodge(jitter.width = 0.75, 
                                                dodge.width = 0.75)) +
     geom_hline(aes(yintercept = 0), color = "darkgrey", lty = 2) +
     
     facet_wrap(vars(module)) +
     theme(legend.position = "none") +
     xlab("") +
     ggtitle("Modules 12, 13, 14")
    
# 14   
   p_14a <- df_kme %>% 
     filter(module == 14 & !is.na(anti_anxiety)) %>% 
     mutate(anti_anxiety = ifelse(anti_anxiety == 1, "anti-anxiety use", "no anti-anxiety use")) %>% 
     ggplot(aes(x = anti_anxiety, y = kme)) +
     geom_boxplot(outlier.shape = NA) +
     geom_point(position = position_jitterdodge(jitter.width = 0.5, 
                                                dodge.width = 0.75),
                aes(shape = anti_anxiety)) +
     scale_shape_manual(values = c(1, 4)) +
     geom_hline(aes(yintercept = 0), color = "darkgrey", lty = 2) +
     labs(x = "") +
     ggtitle("Module 14") +
     theme(legend.position = "none")
   
   p_14b <- df_kme %>% 
     filter(module == 14 & !is.na(anti_anxiety)) %>% 
     mutate(anti_anxiety = ifelse(anti_anxiety == 1, "anti-anxiety use", "no anti-anxiety use")) %>% 
     ggplot(aes(x = anti_anxiety, y = kme, color = dx)) +
     geom_boxplot(outlier.shape = NA) +
     geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                                dodge.width = 0.75),
                aes(shape = anti_anxiety)) + 
     scale_shape_manual(values = c(1, 4)) +
     geom_hline(aes(yintercept = 0), color = "darkgrey", lty = 2) +
     labs(x = "") +
     theme(legend.position = "none")
   
   p12_13_14 / (p_14a + p_14b)
   
   