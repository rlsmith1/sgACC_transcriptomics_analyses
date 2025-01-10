

########################################################################################

# run GO and cell-type enrichment

########################################################################################


# libraries ---------------------------------------------------------------


  library(tidyverse)
  library(patchwork)
  library(janitor)
  library(RColorBrewer)
  library(tidymodels)
  library(tidytext)
  library(stm)
  library(ggwordcloud)


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
  ctrl_prefix <- "20221209_55samplesCTRL_19kgenes_vsd_qSVA17_resids"
  # ctrl_prefix2 <- "20221209_55samplesCTRL_19kgenes_vsd_qSVA17_resids_sft4minSize30cutHeight0.98"
  

# LOAD REGRESSED DATA
  load(paste0("objects/", prefix, ".RDS")) # df_vsd_regress

# FILTER FOR CONTROL ONLY
  df_control <- df_vsd_regress %>% 
    filter(grepl("control", sample))

# LOAD CONTROL MODULE ASSIGNMENTS
  load(paste0("objects/", ctrl_prefix, "_SFT4_MODS.RDS")) # df_modules

# IDENTIFY MODULE SET OF INTEREST  
  df_modules_filt <- df_modules %>% 
    filter(min_size == 50 & cut_height == 0.99) %>% 
    dplyr::select(ensembl_gene_id, module) %>% 
    mutate(module = factor(module, levels = 0:max(module)))
  
  
  
# module functional enrichment --------------------------------------------


# SOURCE FUNCTION
  source("functions/f_top_GO_modules.R")
  load("objects/df_go_full_definitions.Rdata") # full go definitions

# PULL ALL WGCNA GENES AFTER FILTERING FOR BACKGROUND
  l_gene_universe <- df_modules$ensembl_gene_id %>% unique

# CREATE GO ANNOTATION

  # get gene info from ensembl database
  hsapiens_ensembl <- useEnsembl(biomart = "genes",
                                 dataset = "hsapiens_gene_ensembl")
  
  # map each gene to GO term
  m_go_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "namespace_1003"),
                    filters = "ensembl_gene_id",
                    values = l_gene_universe,
                    mart = hsapiens_ensembl)
  
  # build annotation list
  l_gene_2_GO <- m_go_ids[,c(1,3)] %>% unstack


# RUN ENRICHMENT FUNCTION ON ALL MODULES
  doParallel::registerDoParallel()
  
  n_mods <- df_modules_filt$module %>% levels %>% .[length(.)] %>% as.numeric
  df_mods_go <- 0:n_mods %>% 
    map_dfr( ~ f_top_GO_modules(l_gene_module = df_modules_filt %>%
                                  filter(module == .x) %>%
                                  pull(ensembl_gene_id),
                                l_gene_universe = l_gene_universe,
                                m_go_ids = m_go_ids,
                                l_gene_2_GO = l_gene_2_GO,
                                go_ontology = "BP") %>% 
               mutate(module = .x, .before = 1)
             
    ) %>% 
    clean_names() %>%
    dplyr::select(-term) %>%
    left_join(df_go_defs, by = "go_id") %>%
    dplyr::select(module, go_id, term, everything()) %>%
    distinct()

# EXPORT TO EXCEL & SAVE OBJECT
  mod_set <- "sft4minSize50cutHeight0.99"
  
  write.csv(df_mods_go,
            file = paste0("outputs/", ctrl_prefix, sep = "_", mod_set, "_GO_RES.csv"))
  
  save(df_mods_go, file = paste0("objects/", ctrl_prefix, sep = "_", mod_set, "_GO_RES.RDS"))


# PLOT
  df_mods_go %>% 
    
    group_by(module) %>% 
    arrange(module, p_adj) %>% 
    dplyr::top_n(n = 5) %>% 
    
    ggplot(aes(x = -log10(p_adj), y = reorder_within(str_wrap(term, width = 30), -p_adj, module))) +
    geom_point(aes(size = significant/annotated, color = -log10(p_adj))) +
    labs(y = "") +
    facet_wrap(~module, scales = "free") +
    scale_y_discrete(labels = function(x) gsub("\\_.*", "", x)) +
    
    geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
    geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
    geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +
    scale_color_gradientn(colors = brewer.pal(9, "Reds")[3:9]) +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 10),
          strip.text = element_text(size = 10),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12)) +
    ggtitle("Control module enrichment")

  
  

# Topic modeling ----------------------------------------------------------
  
  
  data(stop_words)
  
# UNNEST TOKENS FOR GO TERMS WITHIN EACH GROUP
  df_go_tokens <- df_mods_go %>% 
    filter(!is.na(term)) %>% 
    mutate(module = factor(module, levels = 0:n_mods)) %>% 
    group_by(module) %>% 
    unnest_tokens(word, term)
  
# CREATE A SPARSE MATRIX USING ONLY WORDS THAT APPEAR AT LEAST 3 TIMES
  go_sparse <- df_go_tokens %>%
    count(module, word) %>%
    filter(n > 3 & 
             !(word %in% c(stop_words$word, 
                           "positive", "negative", "regulation", "response"))) %>%
    cast_sparse(module, word, n)
  
# FIT TOPIC MODEL
  topic_model <- stm(go_sparse, K = 5, verbose = FALSE)
  
  summary(topic_model)
  
# TOPIC DESCRIPTIONS:
  # lift = frequency divided by frequency in other topics
  # FREX weights words by frequency and exclusivity to the topic
  p_topic <- tidy(topic_model, matrix = "frex") %>%
    group_by(topic) %>%
    slice_head(n = 10) %>%
    mutate(topic = paste0("topic ", topic)) %>% 
    
    ggplot(aes(label = term)) +
    geom_text_wordcloud(size = 4) +
    facet_wrap(vars(topic), nrow = 1) +
    theme(strip.text = element_text(size = 12))
  
# GENE LIST - TOPIC PROBABILITIES
  group_gamma <- tidy(
    topic_model, 
    matrix = "gamma",
    document_names = rownames(go_sparse)
  ) %>% 
    mutate(module = factor(document, levels = 0:n_mods)) %>% 
    dplyr::select(module, topic, gamma) %>% 
    mutate(topic = factor(topic))
  
  p_gamma <- group_gamma %>%
    ggplot(aes(x = topic, y = gamma, color = topic)) +
    geom_text_repel(aes(alpha = gamma, label = module), size = 6, max.overlaps = 50) +
    labs(y = expression(gamma), x = "") +
    ggtitle("Module strength of association with each topic") +
    theme(legend.position = "none")
  
  (p_gamma / p_topic) + plot_layout(heights = c(2, 1))
  
  
  
  
  
# cell-type enrichment ----------------------------------------------------


# LOAD CELL-TYPE DATA  
  load("objects/cell_type_data.RDS")

# CREATE MODULE-CELL TYPE DF
  df_modules_cell_types <- df_modules_filt %>%
    left_join(df_cell_type) %>%
    filter(!is.na(type))

# HYPERGEOMETRIC OVERLAP

  source("functions/cell_type_hypergeometric_overlap.R")
  
  # create df of all combinations of cell-type and modules
  df_combos <- expand_grid(
    modules = levels(df_modules_cell_types$module),
    type = unique(df_modules_cell_types$type)
  ) %>% 
    filter(grepl("lake", type) & !(grepl("synapse", type)))

  # find intersect of our data and the cell-type expression dataset
  gene_universe_intersect <- inner_join(df_control, df_cell_type) %>% 
    pull(ensembl_gene_id) %>% 
    unique()
  
  # run a hypergeometric function for all combinations of cell-type and modules
  df_mods_celltype_hyper <- map2_dfr(.x = df_combos$modules,
                                     .y = df_combos$type,
                                     .f = ~ f_cell_type_overlap(df_modules_cell_types, 
                                                                module_no = .x, 
                                                                cell_type = .y
                                     )
  ) %>% 
    mutate(q_val = p.adjust(p_val, method = "BH"))


  df_mods_celltype_hyper %>% filter(q_val < 0.05)

  # save results
  save(df_mods_celltype_hyper, 
       file = paste0("objects/", ctrl_prefix, sep = "_", mod_set, "_CELL_TYPE_RES.RDS"))


# PLOT RESULTS

  df_mods_celltype_hyper %>%
    
    # format for plotting
    mutate(label = ifelse(q_val < 0.05, formatC(q_val, format = "e", digits = 2), "")) %>%
    filter(module != 0) %>%
    mutate(module = factor(module, levels = 0:max(as.numeric(.$module)))) %>%
    mutate(q_val = ifelse(q_val < 0.05, q_val, NA)) %>%
    mutate(type = str_remove(type, "_lake")) %>% 
    
    # plot
    ggplot(aes(x = type, y = module,
               fill = -log10(q_val),
               label = label)) +
    geom_tile(aes(width = 0.95, height = 0.95)) +
    scale_fill_gradient(low = "#FFCCCB", high = "red", na.value = "white") +
    xlab("") +
    ylab("") +
    ggtitle("Control module cell type enrichment; sft = 4, minSize = 50, cutHeight = 0.99") +
    geom_text(size = 5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))  




# developmental expression trajectories -----------------------------------


# LOAD EXPRESSION DATA FROM PSYCHENCODE
  df_psychencode_metadata <- read_table("data/eva/Psychencode_brainIDfiles.txt") %>% 
    clean_names()
  
  df_psychencode_expr <- read_table("data/eva/psychencode_scaledlogtransformed_genexpr.txt")
  
# IDENTIFY SAMPLES FROM THE REGION OF INTEREST (NEOCORTEX)  
  neocortex_samples <- df_psychencode_metadata %>% 
    filter(region_superbroad == "Neocortex") %>% 
    pull(id)
  
# FILTER EXPRESSION DATA FOR SAMPLES AND GENES OF INTEREST  
  df_psychencode_expr_filt <- df_psychencode_expr %>% 
    filter(GENE %in% df_modules_filt$ensembl_gene_id) %>% 
    dplyr::select(GENE, all_of(neocortex_samples)) %>% 
    dplyr::rename("ensembl_gene_id" = "GENE")

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
  
# LR MODELING OF SIGNIFICANT PRE- & POST-BIRTH DIFFERENCES
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
    filter(q_val < 0.05) # all significant
  
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
  
  
  
  
  





