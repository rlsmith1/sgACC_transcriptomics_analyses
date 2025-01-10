
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
library(WGCNA)




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

# ADD DRUG CLASSES TO COVARIATES
drugs_of_abuse <- c("smoker", "nicotine", "alcohol", "cannabis", "cocaine", "opioids", "major_stimulants")
mood_meds <- c("anti_anxiety", "mood_stabilizers", "anti_epileptics")
other_psych_meds <- c("anticholinergics", "anti_depressant", "antipsychotics" )
other_drugs <- c("anti_histamines", "non_psychiatric", "other_psychotropic_drug")

df_covariates_all_clean_numeric <- df_covariates_all_clean_numeric %>% 
  
  mutate(
    
    drugs_of_abuse = ifelse(df_covariates_all_clean_numeric %>% 
                              dplyr::select(all_of(drugs_of_abuse)) %>% 
                              rowSums(na.rm = TRUE) != 0, 1, 0),
    mood_meds = ifelse(df_covariates_all_clean_numeric %>% 
                         dplyr::select(all_of(mood_meds)) %>% 
                         rowSums(na.rm = TRUE) != 0, 1, 0),
    other_psych_meds = ifelse(df_covariates_all_clean_numeric %>% 
                                dplyr::select(all_of(other_psych_meds)) %>% 
                                rowSums(na.rm = TRUE) != 0, 1, 0),
    other_drugs = ifelse(df_covariates_all_clean_numeric %>% 
                           dplyr::select(all_of(other_drugs)) %>% 
                           rowSums(na.rm = TRUE) != 0, 1, 0)
    
  ) %>% dplyr::select(-all_of(c(drugs_of_abuse, mood_meds, other_psych_meds, other_drugs)))



# run PCA on each module --------------------------------------------------


# RUN PCA FOR EACH MODULE
df_mods_pca <- df_modules_filt %>% 
  left_join(df_vsd_regress) %>% 
  pivot_wider(id_cols = c(ensembl_gene_id, module), 
              names_from = sample, 
              values_from = resids) %>% 
  group_by(module) %>% 
  nest() %>% 
  
  mutate(pca = map(.x = data,
                   .f = ~ .x %>% 
                     as.data.frame %>%
                     column_to_rownames("ensembl_gene_id") %>% 
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



# eigengenes --------------------------------------------------------------



# GET EIGENGENES
df_kme <- df_mods_pca %>% 
  filter(pc == "pc1") %>% 
  pivot_wider(id_cols = c(module, sample, pc_score), 
              names_from = covariate, 
              values_from = covariate_val) %>% 
  dplyr::rename("kme" = "pc_score") %>% 
  mutate(dx = ifelse(multi_nomial_dx == 1, "control", "case") %>% 
    factor(levels = c("control", "case")),
  .before = 4)

df_kme <- df_kme %>% 
  mutate(group = paste0(dx, 
                        sep = "_", 
                        ifelse(suicide == 1, "suicide", "noSuicide")),
         .before = 4) %>%
  mutate(group = str_replace(group, "NA", "noSuicide")) %>% 
  mutate(group = factor(group, levels = c("control_noSuicide", "case_noSuicide", "case_suicide")))


# PLOT EIGENGENES OF CASE-CONTROL
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

# SUICIDE?
df_kme %>% filter(!is.na(anti_anxiety)) %>% 
  
  ggplot(aes(x = kme, y = group, color = factor(anti_anxiety))) +
  geom_violin() +
  stat_summary(fun = "median", geom = "crossbar", width = 1) +
  geom_point(position = position_jitter(height = 0.25), shape = 1, alpha = 0.7) +
  facet_wrap(vars(module)) +
  labs(y = "") +
  ggtitle("module eigengene distributions by diagnosis in each module") +
  theme(strip.text = element_text(size = 10))



# linear regression of kme ~ case/control ---------------------------------


df_lm_dx <- df_kme %>% 
  group_by(module) %>% 
  nest() %>% 
  mutate(
    
    lm = map(
      
      .x = data,
      .f = ~ lm(kme ~ group + anti_anxiety, data = .x)
      
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
  dplyr::select(module, term, estimate, p_value)




# top down lm -------------------------------------------------------------

technical_covariates <- c("mapped_percent", "five_prime_three_prime_bias", "rin_acsg", 
                          "rn_aextraction_batch", "library_batch", "pmi", "g_cpercent", 
                          "ph")

confounds <- c("age_death", "brain_weight", "height", "weight")


df_kme_regress <- df_kme %>%
  
  dplyr::select(-all_of(technical_covariates)) %>% 
  
  # regress confounds
  nest() %>% 
  mutate(
    
    kme_resids = map(.x = data,
                     .f = ~ lm(kme ~ ., .x %>% dplyr::select(kme, all_of(confounds))
                               ) %>% residuals()
                     )
    
  ) %>% 
  unnest(cols = c(data, kme_resids)) %>% 
  
  # run models on covariates of interest
  dplyr::select(-c(all_of(confounds), kme, sample, bmi, multi_nomial_dx, group, manner)) %>% 
  nest() %>% 
  mutate(
    
    lm = map(data, ~lm(kme_resids ~ ., data = .x)),
    lm_res = map(lm, ~tidy(.x) %>% 
                   clean_names %>% 
                   filter(term != "(Intercept)") %>% 
                   mutate(p_adj = p.adjust(p_value, method = "BH")))
    
  ) %>% 
  unnest(cols = c(lm_res))

df_kme_regress %>% 
  filter(p_value < 0.05) %>% 
  print(n = nrow(.))



# select best model per module using covariates --------------------------------------------


# LIST POTENTIAL COVARIATES 

  covariates <- df_kme_regress %>% 
    filter(p_value < 0.05) %>% 
    pull(term) %>% unique
  
  
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

  doParallel::registerDoParallel()
  # select best recipe for each module based on adj r-squared on test data
  df_mods_best_rec <- map_dfr(
    
    .x = unique(df_kme$module),
    .f = ~ list( {
      
      print(paste0("module ", .x))
      
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

# ROUND 2: ALL COMBOS OF SIGNIFICANT COVARIATES
  sig_covariates <- df_round1_res %>% 
    filter(p_value < 0.05 & !grepl("dx", term) & mod != 0) %>% 
    pull(term) %>% unique

  sig_covariates2 <- paste(sig_covariates, collapse = " + ")

  all_recs <-  map(.x = c(sig_covariates, sig_covariates2),
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
      
      print(paste0("module ", .x))
      
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
  filter(p_value < 0.05) %>% print(n = nrow(.))
write.csv("outputs/20230316_genes_univariate_lm_bestRec.csv",
          row.names = FALSE)













