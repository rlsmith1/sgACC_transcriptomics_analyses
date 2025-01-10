
# !!!!!! RARE !!!!!!

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
  library(ggrepel)


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

# GENE TO TRANSCRIPT MAPPING
load(paste0(base_dir, "objects/df_hsapiens_genome.RDS")) # df_hsapiens_genome

df_gene_to_transcript <- df_hsapiens_genome %>% 
  dplyr::select(transcript_id, transcript_name, gene_id, gene_name) %>% 
  distinct() %>% 
  dplyr::rename("gene_symbol" = "gene_name", "transcript_symbol" = "transcript_name", 
                "ensembl_gene_id" = "gene_id", "ensembl_transcript_id" = "transcript_id"
  ) %>% 
  filter(!is.na(ensembl_transcript_id))


# map transcripts to gene level -------------------------------------------

  df_modules_filt <- df_modules_filt %>% 
    left_join(df_gene_to_transcript) %>% 
    dplyr::select(ensembl_transcript_id, ensembl_gene_id, gene_symbol, module, color) %>% 
    distinct()
    
  df_modules_filt_geneLevel <- df_modules_filt %>% 
    dplyr::select(-ensembl_transcript_id) %>%
    filter(!is.na(ensembl_gene_id)) %>% 
    distinct() 
  
  # 25,000 transcripts map to 18,180 genes
  
# PLOT MODULE SIZES AT TRANSCRIPT & GENE LEVELS
  df_modules_filt %>% 
    dplyr::count(module, color) %>% 
    mutate(type = "transcript level") %>% 
    bind_rows(
      df_modules_filt_geneLevel %>% 
        dplyr::count(module, color) %>% 
        mutate(type = "gene level")
    ) %>% 
    mutate(type = factor(type, levels = c("transcript level", "gene level"))) %>% 
    
    ggplot(aes(x = module, y = n)) +
    geom_col(aes(fill = module)) +
    geom_text(aes(label = n), size = 2, vjust = -0.2) +
    scale_fill_manual(values = c(df_modules_filt %>% 
                                   dplyr::count(module, color) %>% 
                                   pull(color)
    )
    ) +
    facet_wrap(vars(type)) +
    labs(title = "Number of transcripts & genes in each module",
         caption = "25,000 transcripts map to 18,687 unique genes") +
    theme(legend.position = "none")
  
  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/transcripts/RARE/module_sizes", .x), width = 8, height = 4)
  )
  
  

    
# module functional enrichment --------------------------------------------


# SOURCE FUNCTION
  source(paste0(base_dir, "functions/f_top_GO_modules.R"))
  load(paste0(base_dir, "objects/df_go_full_definitions.Rdata")) # full go definitions

# PULL ALL WGCNA GENES AFTER FILTERING FOR BACKGROUND
  gene_universe <- df_modules_filt_geneLevel %>% pull(ensembl_gene_id) %>% unique

# CREATE GO ANNOTATION

  # get gene info from ensembl database
  hsapiens_ensembl <- useEnsembl(biomart = "genes",
                                 dataset = "hsapiens_gene_ensembl")
  
  # map each gene to GO term
  m_go_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id"),
                    filters = "ensembl_gene_id",
                    values = gene_universe,
                    mart = hsapiens_ensembl)
  
  # build annotation list
  l_gene_2_GO <- m_go_ids[,c(1,3)] %>% unstack


# RUN ENRICHMENT FUNCTION ON ALL MODULES
  doParallel::registerDoParallel()
  
  # n_mods <- df_modules_filt$module %>% levels %>% .[length(.)] %>% as.numeric
  df_combos <- df_modules_filt_geneLevel %>% 
    dplyr::select(-ensembl_gene_id, -gene_symbol) %>% 
    distinct() %>% 
    mutate(mod_set = "2_35", .before = 1)
  
  l_mods_go <- map2(
    
    .x =  df_combos %>% pull(mod_set),
    .y = df_combos %>% pull(module),
    
    .f = ~ list(
      
      print(paste0("calculating mod_set ", .x, " module ", .y)),
      
      f_top_GO_modules(l_gene_module = df_modules_filt_geneLevel %>%
                              filter(module == .y) %>%
                              pull(ensembl_gene_id),
                            l_gene_universe = gene_universe,
                            m_go_ids = m_go_ids,
                            l_gene_2_GO = l_gene_2_GO,
                            go_ontology = "BP") %>% 
      mutate(mod_set = .x, .before = 1) %>% 
      mutate(module = .y, .before = 2) 
      
      )
      
  )
    
  df_mods_go <- map(l_mods_go, 2) %>% 
    bind_rows() %>% 
    clean_names() %>%
    dplyr::select(-term) %>%
    left_join(df_go_defs, by = "go_id") %>%
    dplyr::select(mod_set, module, go_id, term, everything()) %>%
    distinct()
  
# EXPORT TO EXCEL & SAVE OBJECT
  write.csv(df_mods_go,
            file = paste0(base_dir, "outputs/tables/", prefix, sep = "_", "SIGNED_SFT", soft_power, "_GO_RES.csv"))
  
  save(df_mods_go, file = paste0(base_dir, "objects/", prefix, sep = "_", "SIGNED_SFT", soft_power, "_GO_RES.RDS"))
  
# PLOT
  df_mods_go %>% 
    filter(!is.na(term)) %>% 
    dplyr::group_by(module) %>% 
    arrange(module, -weight_fisher) %>% 
    mutate(row = row_number()) %>% 
    top_n(n = 5, wt = row) %>% 
    #slice_max(order_by = -weight_fisher, n = 5) %>% 
    
    # plot
    ggplot(aes(x = -log10(p_adj), 
               y = reorder_within(str_wrap(term, width = 35), -log10(p_adj), module))) +
    geom_point(aes(size = significant/annotated, fill = -log10(p_adj + 0.0001)),
               shape = 21, color = "midnightblue") +
    geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
    scale_size_continuous(range = c(2, 5)) +
    scale_fill_gradient(low = "#BFD3E6", high = "midnightblue", na.value = "midnightblue",
                        guide = "none") +
    scale_y_discrete(labels = function(x) {gsub("\\_.*$", "", x)}) +
    facet_wrap(vars(module), scales = "free_y", nrow = 4) +
    labs(y = "", x = "log10(FDR)",
         title = "Gene Ontology enrichment by module") +
    guides(size = guide_legend(nrow = 2, title = "significant/\nannotated")) +
    theme(legend.position = "none")
  
  # save
  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/transcripts/RARE/sft", soft_power, "_GO_RES", .x),
                  width = 18, height = 10)
  )
  
  


  
# Topic modeling ----------------------------------------------------------

  data(stop_words)
  
# UNNEST TOKENS FOR GO TERMS WITHIN EACH GROUP
  df_go_tokens <- df_mods_go %>% 
    filter(p_adj < 0.05 & module != 0) %>% 
    filter(!is.na(term)) %>% 
    group_by(module) %>% 
    unnest_tokens(word, term)
  
# CREATE A SPARSE MATRIX USING ONLY WORDS THAT APPEAR AT LEAST 3 TIMES
  vague_terms <- c("positive", "negative", "regulation", "response", "protein",
                   "pathway", "process", "involved", "signaling")
  go_sparse <- df_go_tokens %>%
    dplyr::count(module, word) %>%
    mutate(module = as.numeric(as.character(module))) %>% 
    anti_join(stop_words) %>% 
    filter(n > 2 & !(word %in% vague_terms)) %>%
    cast_sparse(module, word, n)
  
  dim(go_sparse)
  
# FIT TOPIC MODEL
  set.seed(456)
  topic_model <- stm(go_sparse, K = 5, verbose = FALSE)
  
# TOPIC DESCRIPTIONS:
  # lift = frequency divided by frequency in other topics
  # FREX weights words by frequency and exclusivity to the topic
  p_topic <- tidy(topic_model, matrix = "beta") %>%
    group_by(topic) %>%
    arrange(topic, -beta) %>% 
    top_n(n = 10, wt = beta) %>%
    mutate(topic = paste0("topic ", topic)) %>% 
    filter(!str_detect(term, "cell")) %>% 
    
    ggplot(aes(x = beta, y = reorder_within(term, beta, topic), fill = beta)) +
    geom_col(color = "midnightblue") +
    scale_y_reordered() +
    facet_wrap(vars(topic), nrow = 1, scales = "free_y") +
    scale_fill_gradient(low = "white", high = "midnightblue", guide = "none", limits = c(0, 0.1), na.value = "midnightblue") +
    labs(y = "", x = "Word strength of association with each topic") +
    theme_classic()
  
# GENE LIST - TOPIC PROBABILITIES
  df_col <- df_modules_filt %>% 
    dplyr::select(module, color) %>% 
    distinct()
  colors <- df_col %>% pull(color)
  names(colors) <- df_col %>% pull(module)
  
  group_gamma <- tidy(
    topic_model, 
    matrix = "gamma",
    document_names = rownames(go_sparse)
  ) %>% 
    mutate(module = factor(document)) %>% 
    dplyr::select(module, topic, gamma) %>% 
    mutate(topic = factor(topic))
  
  p_gamma <- group_gamma %>%
    filter(gamma > 0.01 ) %>% 
    ggplot(aes(x = topic, y = gamma, color = module)) +
    geom_point(aes(alpha = gamma)) +
    geom_text_repel(aes(alpha = gamma, label = module), size = 4, 
                    max.overlaps = 50) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), lty = 2, color = "black") +
    labs(y = "Module strength of association \nwith each topic", x = "") +
    theme(legend.position = "none",
          axis.text.x = element_blank())
  
  p_topic / p_gamma + plot_layout(heights = c(1, 2)) + plot_annotation(title = "Topic modeling of GO results")

# save
  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/transcripts/RARE/topic_modeling", .x),
                  width = 10, height = 5)
  )
  
    
# cell-type enrichment ----------------------------------------------------


# LOAD CELL-TYPE DATA  
  load(paste0(base_dir, "objects/cell_type_data.RDS"))

# CREATE MODULE-CELL TYPE DF
  df_modules_cell_types <- df_modules_filt_geneLevel %>%
    left_join(df_cell_type) %>%
    filter(!is.na(type))

# HYPERGEOMETRIC OVERLAP

  source(paste0(base_dir, "functions/cell_type_hypergeometric_overlap.R"))
  
  # create df of all combinations of cell-type and modules
  df_combos <- expand_grid(
    modules = unique(df_modules_cell_types$module),
    type = unique(df_modules_cell_types$type)
  ) %>% 
    filter(grepl("lake", type) & !(grepl("synapse", type)))

  # find intersect of our data and the cell-type expression dataset
  gene_universe_intersect <- inner_join(df_modules_filt_geneLevel, df_cell_type) %>% 
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
    mutate(p_adj = p.adjust(p_val, method = "BH")) %>% 
    mutate(type = str_remove(type, "_lake")) %>% 
    mutate(type = case_when(
      type == "astro" ~ "Astrocyte",
      type == "endo" ~ "Endothelial",
      type == "micro" ~ "Microglia",
      type == "neuro-ex" ~ "Excitatory N",
      type == "neuro-in" ~ "Inhibitory N",
      type == "oligo" ~ "Oligodendrocyte",
      type == "opc" ~ "OPC"
    )
    )
  
  df_mods_celltype_hyper %>% filter(p_adj < 0.05)
  
  # save results
  save(df_mods_celltype_hyper, 
       file = paste0("objects/", prefix, sep = "_", "SIGNED_SFT", soft_power, "_CELL_TYPE_RES.RDS"))



# PLOT RESULTS
  df_mods_celltype_hyper %>%
    filter(module != 0) %>% 
    mutate(color = factor((ifelse(p_val < 0.05, module, NA) %>% as.numeric) - 1),
           label = ifelse(p_adj < 0.05, "*", "")) %>%
    
    ggplot(aes(x = type, y = module)) +
    geom_tile(aes(fill = color, alpha = -log10(p_adj), color = color), 
              width = 0.95, height = 0.95, linewidth = 0.5) +
    geom_text(aes(label = label), vjust = 0.75, size = 5) +
    scale_alpha_continuous(range = c(0.5, 1)) +
    scale_fill_manual(values = colors, na.value = "transparent") +
    scale_color_manual(values = colors, na.value = "transparent") +
    labs(title = "Cell-type enrichment of modules", x = "Cell type", y = "Module",
         caption = "* = FDR < 0.05") +
    theme_classic() +
    theme(legend.position = "none")
  
  ## SAVE   
  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/transcripts/RARE/cell_type", .x),
                  width = 6, height = 5)
  )
  
  
  
# developmental trajectories ----------------------------------------------
  
  
  load(paste0(base_dir, "objects/psychencode_development_expr_data.Rdata")) # df_psychencode_metadata, df_psychencode_expr
  
# IDENTIFY SAMPLES FROM THE REGION OF INTEREST (NEOCORTEX)  
  neocortex_samples <- df_psychencode_metadata %>% 
    filter(region_superbroad == "Neocortex") %>% 
    pull(id)
  
# FILTER EXPRESSION DATA FOR SAMPLES AND GENES OF INTEREST  
  df_psychencode_expr_filt <- df_psychencode_expr %>% 
    filter(GENE %in% df_modules_filt_geneLevel$ensembl_gene_id) %>% 
    dplyr::select(GENE, all_of(neocortex_samples)) %>% 
    dplyr::rename("ensembl_gene_id" = "GENE")
  
# JOIN EXPR WITH MODULE DATA
  df_expr_mod <- df_psychencode_expr_filt %>% 
    left_join(df_modules_filt_geneLevel) %>% 
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
  
# RENAME WINDOWS WITH TIME PERIOD
  li_windows <- c("PCW5-9",
                  "PCW12-13",
                  "PCW16-18",
                  "PCW19-22",
                  "PCW35-PY0.3",
                  "PY0.5-2.5",
                  "PY2.8-10.7",
                  "PY13-19",
                  "PY21-64")
  names(li_windows) = 1:9  
  
# PLOT MEAN EXPRESSION CHANGES
  
  # # assign color vector for modules
  # df_col <- df_modules_filt %>% 
  #   dplyr::select(module, color) %>% 
  #   distinct()
  # colors <- df_col %>% pull(color)
  # names(colors) <- df_col %>% pull(module)
  
  # plot curves
  df_expr_mod_mean %>% 
    ggplot(aes(x = as.numeric(window), y = mean_expr, color = module, fill = module)) +
    geom_smooth(method = "loess") +
    geom_vline(xintercept = 5, lty = 2) +
    scale_x_continuous(breaks = seq(1, length(li_windows), 1),
                       labels = li_windows) +
    facet_wrap( ~ module, scales = "free_y", nrow = 3) +
    
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    
    labs(x = "Developmental window", y = "Mean expression across neocortex samples",
         title = "Average gene expression per module across development") +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 9, angle = 45, hjust = 1))
  
  ## SAVE   
  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/transcripts/RARE/developmental_expression_trajectory", .x),
                  width = 12, height = 6)
  )
  

# put plots together ------------------------------------------------------
  
  theme_set(theme_bw() +
              theme(plot.title = element_text(size = 12),
                    axis.title = element_text(size = 12),
                    axis.text = element_text(size = 10),
                    strip.text = element_text(size = 10)))
  
  # module sizes, cell-type enrichment, topic modeling
  
  ((p_module_sizes + p_cell_type) + plot_layout(widths = c(1.5, 1))) /
    
    ((p_gamma / p_topic) + plot_layout(heights = c(2, 1)))

  