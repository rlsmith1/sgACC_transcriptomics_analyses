

########################################################################################

# GENE LEVEL: run GO and cell-type enrichment

########################################################################################


# setup -------------------------------------------------------------------


## LIBRARIES
library(tidyverse)
library(patchwork)
library(janitor)
library(RColorBrewer)
library(tidymodels)
library(tidytext)
library(stm)
library(ggwordcloud)
library(ggrepel)
library(readxl)
library(writexl)
library(clusterProfiler)

## SET DIRECTORIES
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "08Mar2024_GENES_qSVAgeSexRaceGC"

## IDENTIFY MODULE SET OF INTEREST
soft_power <- 3
minimum_size <- 40
tree_cut_height <- 0.98

## LOAD DATA OBJECTS 
load(paste0(base_dir, "objects/hsapiens_genome_v110.RDS")) # df_hsapiens_genome
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

## FILTER FOR MODULE SET OF INTEREST  
  df_modules_filt <- df_modules %>% 
      filter(min_size == minimum_size & cut_height == tree_cut_height) %>% 
      unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
      arrange(mod_set, module) %>% 
      
      # rename modules to specify gene-level
      mutate(module = paste0("geneM", module) %>% factor(levels = paste0("geneM", levels(module)))
      )
  
 ## ENSEMBL ID TO GENE SYMBOL MAPPINGS
  df_ensembl_to_symbol <- df_hsapiens_genome %>% 
      dplyr::select(gene_id, gene_name) %>% 
      distinct() %>% 
      dplyr::rename("gene_symbol" = "gene_name", "ensembl_gene_id" = "gene_id")
  
  
## PLOT THEME
  theme_set(theme_bw() +
                theme(plot.title = element_text(size = 12),
                      axis.title = element_text(size = 12),
                      axis.text = element_text(size = 10),
                      strip.text = element_text(size = 10),
                      legend.title = element_text(size = 10),
                      legend.text = element_text(size = 10)
                )
  )
  module_colors <- df_modules_filt %>% dplyr::select(module, color) %>% distinct %>% deframe
  
  
  
# Fig 2A | module sizes ------------------------------------------------------------

  df_modules_filt %>% 
    dplyr::count(module, color) %>% 
    ggplot(aes(x = module, y = n)) +
    geom_col(aes(fill = module), color = "black", linewidth = 0.25) +
    geom_text(aes(label = n), size = 2, vjust = -0.2) +
    scale_fill_manual(values = module_colors) +
    labs(title = "Module sizes",
         x = "", y = "Number of genes") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 0.9, size = 7)
          )

# save to project dir
  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/genes/2A.module_sizes", .x),
                  width = 5, height = 4)
  )

# Fig S5 | module functional enrichment --------------------------------------------

  
# RUN GO OVER-REPRESENTATION ANALYSIS ON EACH MODULE
gene_universe <- df_modules_filt %>% pull(ensembl_gene_id) %>% unique 
length(gene_universe) # n = 18677
modules <- df_modules_filt %>% pull(module) %>% unique
  
doParallel::registerDoParallel()
df_mods_go <- map_dfr(
  
  .x = modules,
  .f = ~ enrichGO(gene = df_modules_filt %>% filter(module == .x) %>% pull(ensembl_gene_id),
                  OrgDb = "org.Hs.eg.db",
                  universe = gene_universe,
                  keyType = "ENSEMBL",
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  readable = TRUE
  ) %>% 
    as_tibble() %>% 
    mutate(module = .x, .before = 1)
  
)

save(df_mods_go, file = paste0(base_dir, "objects/", prefix, sep = "_", "SIGNED_SFT", soft_power, "_GO_RES.RDS"))
  
# PLOT RESULTS
load(paste0(base_dir, "objects/", prefix, sep = "_", "SIGNED_SFT", soft_power, "_GO_RES.RDS"))
df_mods_go %>% 
    clean_names %>% 
    filter(!is.na(description) & module != 0) %>% 
    
    # take top 5 paths by p-value per module to plot
    dplyr::group_by(module) %>% 
    arrange(module, -pvalue) %>% 
    mutate(row = row_number(),
           gene_ratio = parse(text = gene_ratio) %>% eval
    ) %>% 
    top_n(n = 5, wt = row) %>% 
    #slice_max(order_by = -weight_fisher, n = 5) %>% 
    mutate(`survives FDR` = ifelse(p_adjust < 0.05, "yes", "no")) %>% 
    mutate(description = str_wrap(description, width = 35)) %>%
    
    # plot
    ggplot(aes(x = -log10(pvalue), y = reorder_within(description, within = module, by = -pvalue))) +
    #y = reorder_within(str_wrap(description, width = 35), -pvalue, module))) +
    geom_point(aes(size = gene_ratio, fill = ontology, color = `survives FDR`),
               shape = 21, stroke = 1) +
    geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
    scale_size_continuous(range = c(2, 5)) +
    scale_color_manual(values = c("yes" = "black", "no" = "transparent")) +
    scale_y_reordered() +
    facet_wrap(vars(module), scales = "free_y", ncol = 4) +
    labs(y = "", x = "log10(FDR)",
         title = "Over-represented GO terms of gene-level modules"
    ) +
    guides(size = guide_legend(title = "Gene ratio")) +
    theme(legend.position = "bottom",
          legend.title.position = "top")

# save to project dir
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/genes/S5.module_GO", .x),
                width = 13, height = 13)
)

# export module assignments and GO results
write_xlsx(list(
    "A. Module gene assignments" = df_modules_filt %>% left_join(df_ensembl_to_symbol) %>% dplyr::select(module, color, ensembl_gene_id, gene_symbol) %>% as.data.frame,
    "B. Module GO results" = df_mods_go %>% filter(pvalue < 0.05 & module != 0) %>% as.data.frame
),
path = paste0(base_dir, "outputs/tables/for_manuscript/TableS1_GENES_module_GO.xlsx")
)



  
# Fig 2B | topic modeling ----------------------------------------------------------

  data(stop_words)
  
# UNNEST TOKENS FOR GO TERMS WITHIN EACH GROUP
  df_go_tokens <- df_mods_go %>% 
    filter(pvalue < 0.05) %>% 
    filter(!is.na(Description)) %>% 
    group_by(module) %>% 
    unnest_tokens(word, Description)
  
# CREATE A SPARSE MATRIX USING ONLY WORDS THAT APPEAR AT LEAST 3 TIMES
  vague_terms <- c("positive", "negative", "regulation", "response", "protein", "activity",
                   "pathway", "process", "involved", "signaling")
  go_sparse <- df_go_tokens %>%
      dplyr::count(module, word) %>%
      mutate(module = as.numeric(as.character(module %>% str_remove("geneM")))) %>% 
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
    geom_col(color = "#222222") +
    scale_y_reordered() +
    scale_x_continuous(limits = c(0, 0.06), breaks = c(0, 0.05)) +
    facet_wrap(vars(topic), nrow = 1, scales = "free_y") +
    scale_fill_gradient(low = "white", high = "#222222", guide = "none", limits = c(0, 0.05), na.value = "#222222") +
    labs(y = "", x = "Word strength of association with each topic") +
    theme_classic()
  
# GENE LIST - TOPIC PROBABILITIES
  group_gamma <- tidy(
      topic_model, 
      matrix = "gamma",
      document_names = rownames(go_sparse)
  ) %>% 
      mutate(module = factor(document)) %>% 
      dplyr::select(module, topic, gamma) %>% 
      mutate(topic = factor(topic),
             module = paste0("geneM", module) %>% factor(levels = paste0("geneM", levels(module)))
      )
  
  p_gamma <- group_gamma %>%
    filter(gamma > 0.01) %>% 
      
    ggplot(aes(x = topic, y = gamma, color = module)) +
    geom_jitter(aes(alpha = gamma), position = position_jitter(seed = 2, width = 0.15)) +
    geom_text_repel(aes(alpha = gamma, label = module), size = 3, box.padding = 0.2, min.segment.length = 0,
                    max.overlaps = 7, position = position_jitter(seed = 2, width = 0.15)) +
    scale_color_manual(values = module_colors) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), color = "black", linewidth = 0.25) +
    labs(y = "Module strength of association \nwith each topic", x = "") +
    theme(legend.position = "none",
          axis.text.x = element_blank())
  
 p_topic / p_gamma + plot_layout(heights = c(1, 2)) + plot_annotation(title = "Semantic themes across module GO results")
  
 
# save to project dir
 map(
   .x = c(".png", ".pdf"),
   .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/genes/2B.topic_modeling", .x),
                 width = 14, height = 5)
 )
 

  
# FIg 2C | cell-type enrichment ----------------------------------------------------

# LOAD CELL-TYPE DATA  
  load(paste0(base_dir, "objects/cell_type_data.RDS"))

# CREATE MODULE-CELL TYPE DF
  df_modules_cell_types <- df_modules_filt %>%
    
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
  gene_universe_intersect <- intersect(colnames(df_vsd_regress), df_cell_type$ensembl_gene_id) %>% 
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

# PLOT RESULTS
  df_mods_celltype_hyper %>%
      mutate(color = ifelse(p_val < 0.05, paste0(module), NA) %>% factor,
             label = ifelse(p_adj < 0.05, "*", "")
      ) %>%
      
      ggplot(aes(x = type, y = module)) +
      geom_tile(aes(fill = color, alpha = -log10(p_adj), color = color), 
                width = 0.95, height = 0.95, linewidth = 0.5) +
      geom_text(aes(label = label), vjust = 0.75, size = 5) +
      scale_alpha_continuous(range = c(0.5, 1)) +
      scale_fill_manual(values = module_colors, na.value = "transparent") +
      scale_color_manual(values = module_colors, na.value = "transparent") +
      labs(title =  "Cell-type enrichment of modules", 
           x = "Cell type", y = "",
           caption = "* = FDR < 0.05") +
      theme_classic() +
      theme(legend.position = "none")
  
# save to project dir
  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/genes/2C.cell_type", .x),
                  width = 6, height = 5)
  )
  

# Fig 2D | developmental trajectories ----------------------------------------------

  load(paste0(base_dir, "objects/psychencode_development_expr_data.Rdata")) # df_psychencode_metadata, df_psychencode_expr
  
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
 
# DEFINE PLOTTING FUNCTION
  p_plot_traj <- function(df,
                          xlab = "Developmental window",
                          ylab = "Mean expression across neocortex samples") {
    df %>% 
      ggplot(aes(x = as.numeric(window), y = mean_expr, color = module, fill = module)) +
      geom_smooth(method = "loess") +
      geom_vline(xintercept = 5, linewidth = 0.2) +
      scale_x_continuous(breaks = seq(1, length(li_windows), 1),
                         labels = li_windows) +
      facet_wrap( ~ module, scales = "free_y", nrow = 3) +
      scale_fill_manual(values = module_colors) +
      scale_color_manual(values = module_colors) +
      labs(x = paste0(xlab), y = paste0(ylab), title = unique(df$group)) +
      theme(legend.position = "none",
            axis.text.x = element_blank())
  }
  
  
  p_plot_traj(df_expr_mod_mean)
  
  # save to project dir
  map(
      .x = c(".png", ".pdf"),
      .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/genes/2D.developmental_expression", .x),
                    width = 14, height = 5)
  )
  
   
# CLUSTER CURVES FOR PLOTTING

  # group 1: low prenatal expression, increase in expression throughout development (M1, M8, M11, M14, M15, M20)
  # group 2: high prenatal expression, decrease in expression throughout development (M2, M12)
  # group 3: low prenatal expression, peak in expression right after birth, then decrease for the rest of development (M3, M5, M6, M16, M18, M19, M21)
  # group 4: increase in expression prenatally, peak right before birth, decrease immediately after birth, then increase for the rest of development (M4, M9, M13, M17, M22, M23)
  
  # other: increase in expression prenatally, peak right before birth, then decrease postnatally
  
  # less clear patterns: M7, M10
  
  module_groups <- c("Category 1", "Category 2", "Category 3", "Category 4",
                     "Category 3", "Category 3",  "other",  "Category 1",
                     "Category 4",  "other",  "Category 1",  "Category 2",
                     "Category 4",  "Category 1",  "Category 1",  "Category 3",
                     "Category 4",  "Category 3",  "Category 3",  "Category 1",
                     "Category 3",  "Category 4",  "Category 4"
  )  
  names(module_groups) <- 1:23

  df_expr_mod_mean_grouped <- df_expr_mod_mean %>% 
      left_join(enframe(module_groups) %>% 
                    dplyr::rename_all(~c("module", "group")) %>% 
                    mutate(module = factor(module, levels = 1:23))
      ) %>% 
      mutate(group = factor(group, levels = c(paste0("Category ", 1:4), "other"))) %>% 
      arrange(group) %>% 
      mutate(module = factor(module, levels = unique(.$module)))
  
  p_plot_traj(df_expr_mod_mean_grouped %>% filter(module != 0)) +
      labs(title = "Average developmental expression pattern")
  
# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/genes/2D.developmental_expression", .x),
                  width = 10, height = 5)
)



#### PLOT AND PATCH TOGETHER
layout <- c(
  patchwork::area(t = 1, b = 8*(180/8), l = 1, r = 60),
  patchwork::area(t = 1, b = 3*(180/8), l = 61, r = 120),
  patchwork::area(t = 3*(180/8) + 1, b = 8*(180/8), l = 61, r = 120),
  patchwork::area(t = 1, b = 3*(180/8), l = 121, r = 180),
  patchwork::area(t = 4*(180/8) + 1, b = 8*(180/8), l = 121, r = 180)
)

p_fig1d <- (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 1"), xlab = "")) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 2"), xlab = "", ylab = "") + theme(axis.text.x = element_blank())) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 3"), ylab = "")) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 4"), xlab = "", ylab = "") + theme(axis.text.x = element_blank())) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 5"), xlab = "", ylab = "")) +
  
  plot_layout(design = layout) +
  plot_annotation(title =" Average developmental expression pattern")

p_fig1d



  