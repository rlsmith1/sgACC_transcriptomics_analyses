

########################################################################################

# OVERALL TRANSCRIPTS: run GO and cell-type enrichment

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
library(ggalluvial)
library(clusterProfiler)

## SET DIRECTORIES
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC_CVq1"

## PARAMETERS FOR MODULE SELECTION
soft_power <- 2
minimum_size <- 35
tree_cut_height <- 0.988

## LOAD OBJECTS 
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

## IDENTIFY MODULE SET OF INTEREST  
df_modules_filt <- df_modules %>% 
    filter(min_size == minimum_size & cut_height == tree_cut_height) %>% 
    unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module) %>% 
    
    # rename modules to specify gene-level
    mutate(module = paste0("transcriptM", module) %>% factor(levels = paste0("transcriptM", levels(module)))
    )

## GENE TO TRANSCRIPT MAPPING
load(paste0(base_dir, "objects/hsapiens_genome_v110.RDS")) # df_hsapiens_genome
df_gene_to_transcript <- df_hsapiens_genome %>% 
  dplyr::select(transcript_id, transcript_name, gene_id, gene_name) %>% 
  distinct() %>% 
  dplyr::rename("gene_symbol" = "gene_name", "transcript_symbol" = "transcript_name", 
                "ensembl_gene_id" = "gene_id", "ensembl_transcript_id" = "transcript_id"
  ) %>% 
  filter(!is.na(ensembl_transcript_id)) %>% 
  mutate(transcript_symbol = ifelse(is.na(transcript_symbol), ensembl_transcript_id, transcript_symbol),
         gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)
  )
df_ensembl_to_symbol <- df_gene_to_transcript %>% 
    dplyr::select(ensembl_transcript_id, transcript_symbol) %>% 
    distinct()

## PLOT THEME
theme_set(theme_bw() +
              theme(plot.title = element_text(size = 12),
                    axis.title = element_text(size = 12),
                    axis.text = element_text(size = 10),
                    strip.text = element_text(size = 10),
                    legend.title = element_text(size = 10),
                    legend.text = element_text(size = 10)))
module_colors <- df_modules_filt %>% dplyr::select(module, color) %>% distinct %>% deframe


# map transcripts to gene level -------------------------------------------

df_modules_filt <- df_modules_filt %>% 
  left_join(df_gene_to_transcript) %>% 
  dplyr::select(ensembl_transcript_id, ensembl_gene_id, gene_symbol, module, color) %>% 
  distinct()

df_modules_filt_geneLevel <- df_modules_filt %>% 
  dplyr::select(-ensembl_transcript_id) %>%
  filter(!is.na(ensembl_gene_id)) %>% 
  distinct() 

n_transcripts <- df_modules_filt %>% pull(ensembl_transcript_id) %>% unique %>% length # n = 58775
n_genes <- df_modules_filt_geneLevel %>% pull(ensembl_gene_id) %>% unique %>% length # n = 17250

print(paste0(n_transcripts, " transcripts map to ", n_genes, " genes"))


# Fig S7A | PLOT MODULE SIZES AT TRANSCRIPT & GENE LEVELS -------------------------------------------

# COMBINE TRANSCRIPT AND GENE COUNTS FOR EACH MODULE
df_fig5Sa <- df_modules_filt %>% 
  dplyr::count(module, color) %>% 
  mutate(type = "transcript level") %>% 
  bind_rows(
    df_modules_filt_geneLevel %>% 
      dplyr::count(module, color) %>% 
      mutate(type = "gene level")
  ) %>% 
  mutate(type = factor(type, levels = c("transcript level", "gene level")))
  
## PLOT
p_figS5a_transcript <- df_fig5Sa %>% 
    filter(type == "transcript level") %>% 
    ggplot(aes(x = n, y = module)) +
    geom_col(aes(fill = I(color)), color = "black") +
    geom_text(aes(label = n), size = 2.5, vjust = 0.5, hjust = 1) +
    scale_y_discrete(position = "right") +
    scale_x_reverse() +
    labs(title = "Number of transcripts in each module", y = "") +
    theme(axis.text.y = element_text(hjust = -1))

p_figS5a_gene <- df_fig5Sa %>% 
    filter(type == "gene level") %>% 
    ggplot(aes(x = n, y = module)) +
    geom_col(aes(fill = I(color)), color = "black") +
    geom_text(aes(label = n), size = 2.5, vjust = 0.5, hjust = -0.15) +
    labs(y = "", title = "Number of unique corresponding genes in each module",
         caption = paste0(n_transcripts, " transcripts map to ", n_genes, " unique genes")) +
    theme(axis.text.y = element_blank())

p_figS5a_transcript + plot_spacer() + p_figS5a_gene + plot_layout(widths = c(4, -0.3, 4.5), guides = "collect") #+
#plot_annotation(title = "Number of transcripts (and corresponding genes) in each module")

# save to project dir
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/S7A.module_sizes", .x),
                width = 13, height = 5)
)




# Fig S8 | module functional enrichment --------------------------------------------

# RUN GO OVER-REPRESENTATION ANALYSIS ON EACH MODULE
gene_universe <- df_modules_filt_geneLevel %>% pull(ensembl_gene_id) %>% unique
modules <- df_modules_filt_geneLevel %>% pull(module) %>% unique

doParallel::registerDoParallel()
df_mods_go <- map_dfr(
  
  .x = modules,
  .f = ~ enrichGO(gene = df_modules_filt_geneLevel %>% filter(module == .x) %>% pull(ensembl_gene_id),
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

save(df_mods_go, file = paste0(base_dir, "objects/", prefix,
                               "_SIGNED_SFT", soft_power, 
                               "_MIN_SIZE", minimum_size, 
                               "_CUT_HEIGHT", tree_cut_height, "_GO_RES.RDS")
     )

# PLOT RESULTS

load( paste0(base_dir, "objects/", prefix,
             "_SIGNED_SFT", soft_power, 
             "_MIN_SIZE", minimum_size, 
             "_CUT_HEIGHT", tree_cut_height, "_GO_RES.RDS")
)
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
    labs(y = "", x = "log10(p-value)",
         title = "Transcript-level modules GO enrichments") +
    guides(size = guide_legend(title = "Gene ratio")) +
    theme(legend.position = c(0.6, 0.05),
          legend.box = "horizontal",
          legend.direction = "horizontal",
          legend.title.position = "top"
    )

# save to project dir
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/S8.GO_RES", .x),
                width = 13.5, height = 13)
)

# save for publication   
# export module assignments and GO results
write_xlsx(list(
    "A. Module gene assignments" = df_modules_filt %>% left_join(df_ensembl_to_symbol) %>% 
        dplyr::select(module, color, ensembl_transcript_id, transcript_symbol, ensembl_gene_id, gene_symbol) %>% as.data.frame,
    "B. Module GO results" = df_mods_go %>% filter(pvalue < 0.05 & module != 0) %>% as.data.frame
),
path = paste0(base_dir, "outputs/tables/for_manuscript/TableS5_TRANSCRIPTS_module_GO.xlsx")
)



# Fig S7B | topic modeling ----------------------------------------------------------

## LOAD DATA
load( paste0(base_dir, "objects/", prefix,
             "_SIGNED_SFT", soft_power, 
             "_MIN_SIZE", minimum_size, 
             "_CUT_HEIGHT", tree_cut_height, "_GO_RES.RDS")
)
data(stop_words)

# UNNEST TOKENS FOR GO TERMS WITHIN EACH GROUP
df_go_tokens <- df_mods_go %>% 
  filter(p.adjust < 0.05) %>% 
  group_by(module) %>% 
  unnest_tokens(word, Description)

# CREATE A SPARSE MATRIX USING ONLY WORDS THAT APPEAR AT LEAST 3 TIMES
vague_terms <- c("positive", "negative", "regulation", "response", "protein",
                 "pathway", "process", "involved", "signaling")
go_sparse <- df_go_tokens %>%
  dplyr::count(module, word) %>%
  mutate(module = as.numeric(as.character(module %>% str_remove("transcriptM")))) %>% 
  anti_join(stop_words) %>% 
  filter(n > 2 & !(word %in% vague_terms)) %>%
  cast_sparse(module, word, n)

dim(go_sparse)

# FIT TOPIC MODEL
n_topics <- 5
set.seed(456)
topic_model <- stm(go_sparse, K = n_topics, verbose = FALSE)

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
  scale_x_continuous(limits = c(0, 0.11), breaks = c(0, 0.05, 0.1)) +
  facet_wrap(vars(topic), nrow = 1, scales = "free_y") +
  scale_fill_gradient(low = "white", high = "#222222", guide = "none", limits = c(0, 0.1), na.value = "#222222") +
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
           module = paste0("transcriptM", module) %>% factor(levels = paste0("transcriptM", levels(module)))
    )

p_gamma <- group_gamma %>%
    filter(gamma > 0.01 ) %>% 
    ggplot(aes(x = topic, y = gamma, color = module)) +
    geom_jitter(aes(alpha = gamma), position = position_jitter(seed = 2, width = 0.15)) +
    geom_text_repel(aes(alpha = gamma, label = module), size = 2.5, box.padding = 0.2, min.segment.length = 0,
                    max.overlaps = 7, position = position_jitter(seed = 2, width = 0.15)) +
    geom_vline(xintercept = seq(1, n_topics - 1) + 0.5, lty = 2, color = "black") +
    scale_color_manual(values = module_colors) +
    labs(y = "Module strength of association \nwith each topic", x = "") +
    theme(legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

p_topic / p_gamma + plot_layout(heights = c(1, 2)) + plot_annotation(title = "Semantic themes across transcript-level module GO results")

# save to project dir
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/S7B.module_topic_modeling", .x),
                width = 9, height = 4.5)
)



# Fig S7C | cell-type enrichment ----------------------------------------------------


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


# PLOT RESULTS
df_mods_celltype_hyper %>%
  filter(module != 0) %>% 
  mutate(color = ifelse(p_val < 0.05, paste0(module), NA) %>% factor,
         label = ifelse(p_adj < 0.05, "*", "")) %>%
  
  ggplot(aes(x = type, y = module)) +
  geom_tile(aes(fill = color, alpha = -log10(p_adj), color = color), 
            width = 0.95, height = 0.95, linewidth = 0.5) +
  geom_text(aes(label = label), vjust = 0.75, size = 5) +
  scale_alpha_continuous(range = c(0.5, 1)) +
  scale_fill_manual(values = module_colors, na.value = "transparent") +
  scale_color_manual(values = module_colors, na.value = "transparent") +
  labs(title = "Cell-type enrichment", 
       x = "Cell type", y = "",
       caption = "* = FDR < 0.05") +
  theme_classic() +
  theme(legend.position = "none")

# save to project dir
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/S7C.module_cell_type", .x),
                width = 6, height = 5.5)
)



# Fig S7D | developmental trajectories ----------------------------------------------


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


# DEFINE PLOTTING FUNCTION
p_plot_traj <- function(df,
                        xlab = "Developmental window",
                        ylab = "Mean expression across neocortex samples") {
    df %>% 
        ggplot(aes(x = as.numeric(window), y = mean_expr, color = module, fill = module)) +
        geom_smooth(method = "loess") +
        geom_vline(xintercept = 5, lty = 2) +
        scale_x_continuous(breaks = seq(1, length(li_windows), 1),
                           labels = li_windows) +
        facet_wrap( ~ module, scales = "free_y", nrow = 3) +
        scale_fill_manual(values = module_colors) +
        scale_color_manual(values = module_colors) +
        labs(x = paste0(xlab), y = paste0(ylab), title = unique(df$group)) +
        theme(legend.position = "none",
              axis.text.x = element_blank()
        )
}

# PLOT
p_plot_traj(df_expr_mod_mean) + labs(title = "Average gene expression across development")

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/S7D.module_developmental_trajectories", .x),
                  width = 14, height = 5)
)


# **************************************************************************************************

### CLUSTER CURVES FOR PLOTTING

# group 1: low prenatal expression, increase in expression throughout development (M1, M3, M5, M8, M12, M16, M17, M19, M20, M22)
# group 2: high prenatal expression, decrease in expression throughout development (M23)
# group 3: low prenatal expression, peak in expression right after birth, then decrease for the rest of development (M9, M10, M13, M14)
# group 4: increase in expression prenatally, peak right before birth, then decrease postnatally (M4, M6, M11)
# group 5: increase in expression prenatally, peak right before birth, decrease immediately after birth, then increase for the rest of development (M2, M7, M15, M18, M21)

module_groups <- c("Category 1", "Category 5", "Category 1", "Category 4", "Category 1",
                   "Category 4", "Category 5", "Category 1", "Category 3", "Category 3",
                   "Category 4", "Category 1", "Category 3", "Category 3", "Category 5",
                   "Category 1", "Category 1", "Category 5", "Category 1", "Category 1", 
                   "Category 5", "Category 1", "Category 2"
)  
names(module_groups) <- 1:24

df_expr_mod_mean_grouped <- df_expr_mod_mean %>% 
  #group_by(module) %>% 
  #mutate(z_score = scale(mean_expr)[,1]) %>% 
  left_join(enframe(module_groups) %>% dplyr::rename_all(~c("module", "group"))) %>% 
  arrange(group) %>% 
  mutate(module = factor(module, levels = unique(.$module)))


# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/1E.module_developmental_trajectories", .x),
                  width = 9, height = 7)
)


# ************************************************************

# PLOT AND PATCH TOGETHER LANDSCAPE ORIENTATION 
layout <- c(
  patchwork::area(l = 1, r = 180, t = 1, b = 60),
  patchwork::area(l = 1, r = 1*(180/10) + 1, t = 61, b = 120),
  patchwork::area(l = 5.5*(180/10), r = 10*(180/10), t = 61, b = 120),
  patchwork::area(l = 1, r = 3.5*(180/10), t = 121, b = 180),
  patchwork::area(l = 5*(180/10) + 1, r = 10*(180/10), t = 121, b = 180)
)

p_figS5e <- (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 1")) + theme(axis.text.x = element_blank())) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 2")) + theme(axis.text.x = element_blank())) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 3")) + theme(axis.text.x = element_blank())) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 4"))) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 5"))) +
  
  plot_layout(design = layout) +
  plot_annotation(title = bquote(bold("E")~" | Average gene expression across development")) #&
  #xlab("Developmental window") &
  #ylab("Mean expression across neocortex samples")

p_figS5e

# save for publication   
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0("~/Documents/PhD/manuscripts/BiolPsych2024/figures/supplement/FigS5E", .x),
                width = 15, height = 7.5)
)


# PLOT AND PATCH TOGETHER PORTRAIT ORIENTATION 
layout <- c(
  patchwork::area(t = 1, b = 180, l = 1, r = 60),
  patchwork::area(t = 1, b = 1*(180/10) + 1, l = 61, r = 120),
  patchwork::area(t = 4*(180/10), b = 10*(180/10), l = 61, r = 120),
  patchwork::area(t = 1, b = 6*(180/10), l = 121, r = 180),
  patchwork::area(t = 6*(180/10) + 1, b = 10*(180/10), l = 121, r = 180)
)

p_fig1d <- (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 1"), xlab = "")) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 2"), xlab = "", ylab = "") + theme(axis.text.x = element_blank())) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 3"), ylab = "")) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 4"), xlab = "", ylab = "") + theme(axis.text.x = element_blank())) +
  (p_plot_traj(df_expr_mod_mean %>% filter(group == "Category 5"), xlab = "", ylab = "")) +
  
  plot_layout(design = layout) +
  plot_annotation(title = "Average gene expression per module across development")

p_fig1d

# save for publication   
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0("~/Documents/PhD/manuscripts/BiolPsych2024/figures/Fig4D", .x),
                width = 6, height = 10)
)
