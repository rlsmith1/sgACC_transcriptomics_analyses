
########################################################################################

# identify windows of development where gene modules are over-/under-expressed

########################################################################################


# libraries ---------------------------------------------------------------

  library(tidyverse)
  library(janitor)
  library(tidymodels)
  library(RColorBrewer)


# set theme for plots -----------------------------------------------------

  theme_set(theme_bw() +
              theme(plot.title = element_text(size = 18),
                    axis.title = element_text(size = 15),
                    axis.text = element_text(size = 12),
                    strip.text = element_text(size = 15),
                    legend.title = element_text(size = 15),
                    legend.text = element_text(size = 12)))
  


# data --------------------------------------------------------------------

prefix <- "20221209_185samples_19kgenes_vsd_qSVA17_resids"


# IDENTIFY MODULE SET OF INTEREST
  load(paste0("objects/", prefix, "_SFT3_MODS.RDS")) # df_modules
  df_modules_filt <- df_modules %>% 
    filter(min_size == 50 & cut_height == 0.99) %>% 
    dplyr::select(ensembl_gene_id, module) %>% 
    mutate(module = factor(module, levels = 0:max(module)))

# PsychENCODE DATA
  df_psychencode_metadata <- read_table("data/eva/Psychencode_brainIDfiles.txt") %>% 
    clean_names()
  
  df_psychencode_expr <- read_table("data/eva/psychencode_scaledlogtransformed_genexpr.txt")
  
  # save as object for future use
  save(df_psychencode_metadata, df_psychencode_expr,
       file = "objects/psychencode_development_expr_data.Rdata") #df_psychencode_metadata, df_psychencode_expr
  

# format PsychENCODE data -------------------------------------------------


# IDENTIFY SAMPLES FROM THE REGION OF INTEREST (NEOCORTEX)  
  neocortex_samples <- df_psychencode_metadata %>% 
    filter(region_superbroad == "Neocortex") %>% 
    pull(id)
  
# FILTER EXPRESSION DATA FOR SAMPLES AND GENES OF INTEREST  
  df_psychencode_expr_filt <- df_psychencode_expr %>% 
    filter(GENE %in% df_modules_filt$ensembl_gene_id) %>% 
    dplyr::select(GENE, all_of(neocortex_samples)) %>% 
    dplyr::rename("ensembl_gene_id" = "GENE")
  
  

# calculate mean expression of each module for each time window -----------


# JOIN EXPR WITH MODULE DATA
  df_expr_mod <- df_psychencode_expr_filt %>% 
    left_join(df_modules_filt) %>% 
    dplyr::select(ensembl_gene_id, module, everything()) %>% 
    pivot_longer(starts_with("HS"), names_to = "id", values_to = "expression") %>% 
    left_join(df_psychencode_metadata %>% 
                dplyr::select(id, window) %>% 
                filter(id %in% neocortex_samples)) %>% 
    dplyr::select(ensembl_gene_id, module, id, window, expression) %>% 
    dplyr::rename("sample" = "id") %>% 
    mutate(window = as.factor(window)) 

# FIND MEAN EXPRESSION VALUE OF EACH SAMPLE FOR EACH MODULE  
  df_expr_mod_mean <- df_expr_mod %>% 
    group_by(module, sample, window) %>%
    summarise(mean_expr = mean(expression))
  
# MODEL ASSOCIATIONS
  
  # linear
  df_lm_mod_dev_traj <- levels(df_expr_mod_mean$module) %>% 
    map_dfr( ~ lm(mean_expr ~ window, 
                  data = df_expr_mod_mean %>% filter(module == .x)) %>% 
               tidy %>% 
               mutate(module = .x, .before = 1)) %>% 
    clean_names %>% 
    filter(term != "(Intercept)") %>% 
    mutate(q_val = p.adjust(p_value, method = "BH"))
  
  df_lm_mod_dev_traj %>% 
    filter(q_val < 0.05) %>% 
    print(n = nrow(.))
  
  # logistic
  df_lr_mod_dev_traj <- levels(df_expr_mod_mean$module) %>% 
    map_dfr( ~ glm(period ~ mean_expr,
                   family = "binomial", 
                  data = df_expr_mod_mean %>% 
                    mutate(period = ifelse(as.numeric(as.character(window)) < 5, 0, 1)) %>% 
                    filter(module == .x)) %>% 
               tidy %>% 
               mutate(module = .x, .before = 1)) %>% 
    clean_names %>% 
    filter(term != "(Intercept)") %>% 
    mutate(q_val = p.adjust(p_value, method = "BH"))
  
  df_lr_mod_dev_traj_sig <- df_lr_mod_dev_traj %>% 
    filter(q_val < 0.05)
  
# PLOT MEAN EXPRESSION CHANGES
  
  # assign color vector for modules
  set.seed(789)
  qual_col_pals <- brewer.pal.info %>% filter(category == "qual" & colorblind == TRUE)
  col_vector <-  sample(
    unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))),
    (length(df_x_coefs$module %>% unique) - 1)
  )
  
  # plot curves
  df_expr_mod_mean %>% 
    filter(module %in% df_lr_mod_dev_traj_sig$module & module != 0) %>% 
    
    ggplot(aes(x = as.numeric(window), y = mean_expr, color = module, fill = module)) +
    #geom_point() +
    geom_smooth(method = "loess") +
    geom_vline(xintercept = 5, lty = 2) +
    scale_x_continuous(breaks = seq(1,14,1)) +
    facet_wrap( ~ module, scales = "free") +
    
    scale_fill_manual(values = col_vector) +
    scale_color_manual(values = col_vector) +
    
    labs(x = "Developmental window", y = "Mean expression across neocortex samples") +
    ggtitle("Average gene expression per module across development")
  

  
  
    
  
  