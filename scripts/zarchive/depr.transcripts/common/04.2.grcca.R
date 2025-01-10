
########################################################################################

# Common transcript level GRCCA results from Agoston's toolkit

########################################################################################

# libraries ---------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(janitor)
library(RColorBrewer)
library(biomaRt)
library(tidymodels)  
library(ggrepel)
library(raveio)
library(ggpubr)
#library(ggchicklet)


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
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "20230314_185samples_26kCOMMONtranscripts_vst_qSVA123567_MP_RNAbatch_Race_resids"

# LOAD OBJECTS  
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress_common
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

# IDENTIFY MODULE SET OF INTEREST  
df_modules_filt <- df_modules %>% 
  filter(min_size == 50 & cut_height == 0.99) %>% 
  unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
  arrange(mod_set, module)

# COVARIATES
df_covariates <- df_vsd_regress_common %>% 
  ungroup %>% 
  dplyr::select(-ensembl_transcript_id, -vsd, -resids) %>% 
  distinct()

# GENE SYMBOL TO ENSEMBL GENE ID MAPPING
load(paste0(base_dir, "objects/df_hsapiens_genome.RDS"))

df_ensembl_to_symbol <- df_hsapiens_genome %>% 
  dplyr::select(gene_id, gene_name, transcript_id, transcript_name) %>% 
  distinct() %>% 
  dplyr::rename("gene_symbol" = "gene_name", "ensembl_gene_id" = "gene_id",
                "transcript_symbol" = "transcript_name", "ensembl_transcript_id" = "transcript_id") %>% 
  filter(!is.na(ensembl_transcript_id))


# define drug classes to go in matrix -------------------------------------

n_samples <- df_covariates %>% pull(sample) %>% unique %>% length

psych_drugs <- c("anti_anxiety", "anti_depressant", "antipsychotics", "mood_stabilizers", 
                 "anti_epileptics", "anti_histamines", "anticholinergics", "other_psychotropic_drug")
rec_drugs <- c("smoker", "nicotine", "cocaine", "alcohol", "opioids", "cannabis", 
               "major_stimulants", "non_psychiatric")

df_drugs <- df_covariates %>% 
  ungroup %>% 
  dplyr::select(sample, all_of(c(psych_drugs, rec_drugs))) %>% 
  distinct() %>% 
  
  # combine drug classes that are essentially the same
  mutate(anti_histamines = ifelse(anti_histamines == 1 | anticholinergics == 1, 1, 0),
         nicotine = ifelse(nicotine == 1 | smoker == 1, 1, 0)
  ) %>% 
  dplyr::select(-anticholinergics, -smoker) %>% 
  
  # remove cocaine because only one person was using (outsized effect)
  dplyr::select(-cocaine)

# count NAs for each drug
prop_na <- (df_drugs %>% 
              dplyr::select(-sample) %>% 
              is.na %>% 
              colSums) / n_samples


# define matrices ---------------------------------------------------------------

# X 
df_grcca_x <- df_vsd_regress_common %>% 
  dplyr::select(ensembl_transcript_id, sample, resids) %>% 
  left_join(df_modules_filt) %>% 
  arrange(module) 

X.mat <- df_grcca_x %>% 
  pivot_wider(id_cols = sample, names_from = ensembl_transcript_id, values_from = resids) %>% 
  as.data.frame %>% 
  column_to_rownames("sample")

# X group vector
x_group <- df_grcca_x %>% 
  dplyr::select(ensembl_transcript_id, module) %>% 
  distinct() %>% 
  pull(module) %>% 
  as.numeric

# Y and C (confound) matrices
YC.mat <- df_covariates %>% 
  
  # combine drug classes that are essentially the same
  mutate(anti_histamines = ifelse(anti_histamines == 1 | anticholinergics == 1, 1, 0),
         nicotine = ifelse(nicotine == 1 | smoker == 1, 1, 0)
  ) %>% 
  dplyr::select(-anticholinergics, -smoker) %>% 
  
  # select drugs that don't have too many NAs
  dplyr::select(sample, brain_weight, age_death, all_of(names(prop_na[prop_na < 0.05]))) %>% 
  
  # if drug use is NA, assume 0
  mutate_if(is.numeric, ~ ifelse(is.na(.x), 0, .x)) %>% 
  
  # create dx columns
  mutate(control = ifelse(grepl("control", sample), 1, 0),
         bipolar = ifelse(grepl("bipolar", sample), 1, 0),
         mdd = ifelse(grepl("mdd", sample), 1, 0),
         schizo = ifelse(grepl("schizo", sample), 1, 0)
  ) %>% 
  dplyr::select(sample, bipolar, mdd, schizo, everything(), -control) %>% 
  column_to_rownames("sample")

Y.mat <- YC.mat %>% dplyr::select(-c(brain_weight, age_death))

C.mat <- YC.mat %>% dplyr::select(c(brain_weight, age_death))

# EXPORT FOR AGOSTON'S TOOLBOX
type <- "COMMON/"
project_dir <- "DrugsLessThan9NA_noCoc_regressBrainweightAgedeath/"
cca_data_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir, "data/")
dir.create(cca_data_dir, recursive = TRUE)

write.table(X.mat, file = paste0(cca_data_dir, "X.txt"))
write.table(Y.mat, file = paste0(cca_data_dir, "Y.txt"))
write.table(C.mat, file = paste0(cca_data_dir, "C.txt"))
write.table(x_group, file = paste0(cca_data_dir, "XGroup.txt"))

# GRCCA labels files
df_labels_x <- tibble(ensembl_transcript_id = colnames(X.mat)) %>% 
  left_join(df_modules_filt) %>% 
  left_join(df_ensembl_to_symbol, by = join_by(ensembl_transcript_id)) %>% 
  mutate(Label = row_number(), .before = 1) %>% 
  dplyr::select(Label, module, ensembl_transcript_id, transcript_symbol) %>% 
  dplyr::rename_all(~c("Label", "Category", "Ensembl_id", "HGNC_id")) %>% 
  mutate(HGNC_id = ifelse(is.na(HGNC_id), Ensembl_id, HGNC_id))

df_labels_y <- tibble(Label = colnames(Y.mat))

write.csv(df_labels_x, paste0(cca_data_dir, "LabelsX.csv"), row.names = FALSE)
write.csv(df_labels_y, paste0(cca_data_dir, "LabelsY.csv"), row.names = FALSE)



# GRCCA results -----------------------------------------------------------

# GRCCA RES  
type <- "COMMON/"
project_dir <- "DrugsLessThan9NA_noCoc_regressBrainweightAgedeath/"
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)
#analysis_dir <- "grcca_holdout1-0.20_subsamp5-0.20"
analysis_dir <- "grcca_permutation_VARx0.1_1_mu0.1_lambda_0.9999"
model_1 <- read_mat(paste0(cca_dir, "framework/", analysis_dir, "/res/level1/model_1.mat"))

# READ RESULTS TABLE TO SELECT SPLIT
df_results <- read_table(paste0(cca_dir, "framework/", analysis_dir, "/res/level1/results_table.txt"))
df_results %>% arrange(pval, -correl)

# READ IN RES
split <- 5

df_labels_x <- read_csv(paste0(cca_dir, "data/LabelsX.csv")) %>% 
  mutate(Category = ifelse(Category < 10, paste0("0", Category), as.character(Category)))

df_x_res_all <- read_csv(paste0(cca_dir, "framework/", analysis_dir, "/res/level1/weightX_split", split, ".csv")) %>% 
  left_join(df_labels_x) %>%
  clean_names %>% 
  dplyr::rename("ensembl_transcript_id" = "ensembl_id",
                "transcript_symbol" = "hgnc_id") %>% 
  left_join(df_modules_filt) %>% 
  dplyr::select(-category, -label)

df_y_res_all <- read_csv(paste0(cca_dir, "framework/", analysis_dir, "/res/level1/weightY_split", split, ".csv")) %>% 
  clean_names %>% 
  arrange(-weight)  %>% 
  mutate(label = case_when(
    label == "mdd" ~ "MDD",
    label == "bipolar" ~ "BD", 
    label == "schizo" ~ "SCZ",
    TRUE ~ label
  )) %>% 
  mutate(label = factor(label, levels = .$label))


# LOAD SUBSET RESULTS
df_x_res_subset <- read_csv(paste0(cca_dir, "framework/", analysis_dir, "/res/level1/weightX_subset_split", split, ".csv")) %>% 
  left_join(df_labels_x) %>%
  clean_names %>% 
  dplyr::rename("ensembl_transcript_id" = "ensembl_id",
                "transcript_symbol" = "hgnc_id") %>% 
  left_join(df_modules_filt) %>% 
  dplyr::select(-category, -label)

df_y_res_subset <- read_csv(paste0(cca_dir, "framework/", analysis_dir, "/res/level1/weightY_subset_split", split, ".csv")) %>% 
  clean_names %>% 
  arrange(-weight)  %>% 
  mutate(label = case_when(
    label == "mdd" ~ "MDD",
    label == "bipolar" ~ "BD", 
    label == "schizo" ~ "SCZ",
    TRUE ~ label
  )) %>% 
  mutate(label = factor(label, levels = .$label))

# Y RES
df_y_res_all %>% 
  left_join(
    df_y_res_subset %>% 
      mutate(significant = ifelse(weight == 0, "", "*")) %>% 
      dplyr::select(-weight),
    by = join_by(label)
  ) %>% 
  # mutate(label = str_replace(label, "major_", "major "),
  #        label = str_replace(label, "anti_", "anti-")) %>% 
  ggplot(aes(x = reorder(label, weight), y = weight)) +
  geom_col(aes(fill = label, alpha = abs(weight)), color = "black") +
  geom_text(aes(label = significant), size = 10, vjust = 0) +
  geom_hline(aes(yintercept = 0), color = "black") +
  scale_fill_manual(values = c(dx_colors, drug_colors)) +
  scale_alpha_continuous(range = c(0.1, 1)) +
  scale_x_discrete(labels = function(x) {str_replace(x, "_", " ")}) +
  ylim(c(-0.1, 0.1)) +
  labs(x = "", title = "Y coefficient weights") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 0.9))

# SAVE for supplement
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/common_transcript_yres", .x),
                width = 3.5, height = 3.5)
)


# ============================================================================


# X & Y MATRICES  
  Y.mat <- read.table(paste0(cca_dir, "data/Y.txt"))
  df_y <- Y.mat %>% 
    rownames_to_column("sample") %>% 
    as_tibble()
  
  X.mat <- read.table(paste0(cca_dir, "data/X.txt"))
  df_x <- X.mat %>% 
    t() %>% 
    as.data.frame %>% 
    rownames_to_column("gene_id") %>% 
    as_tibble()

# X & Y RESULTS
  df_labels_x <- read_csv(paste0(cca_dir, "data/LabelsX.csv")) %>% 
    mutate(Category = ifelse(Category < 10, paste0("0", Category), as.character(Category)))
  
  df_x_res <- read_csv(paste0(cca_dir, "framework/grcca_holdout1-0.20_subsamp5-0.20/res/level1/weightX_alpha0.05_split1.csv")) %>% 
    left_join(df_labels_x) %>%
    clean_names %>% 
    dplyr::rename("ensembl_gene_id" = "ensembl_id",
                  "gene_symbol" = "hgnc_id") %>% 
    left_join(df_modules_filt) %>% 
    dplyr::select(-category, -label)
  
  df_y_res <- read_csv(paste0(cca_dir, "framework/grcca_holdout1-0.20_subsamp5-0.20/res/level1/weightY_alpha0.05_split1.csv")) %>% 
    clean_names %>% 
    mutate(covariate = factor(label, levels = c("bipolar", "mdd", "schizo", 
                                                    "opioid", "anti-anxiety"))) %>% 
    dplyr::select(-label)
  
  
# plot Effect 1 results -----------------------------------------------------------
  
  
# EFFECT 1 X-Y COR
  df_lvs <- tibble(lvx = as.matrix(X.mat) %*% model_1$wX,
                    lvy = as.matrix(Y.mat) %*% model_1$wY)
  
  df_lvs %>% 
    ggplot(aes(x = lvx, y = lvy)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = lm, lty = 2, color = "midnightblue") +
    xlab("X latent variable (x matrix • x weights)") +
    ylab("Y latent variable (y matrix • y weights)") +
    #stat_regline_equation() +
    stat_cor(label.y = -0.02, label.x = -5.2) +
    ggtitle("Gene x-y latent variable correlation")
  
  
# LEVEL 1 Y   
  df_y_res_overall <- tibble(covariate = colnames(Y.mat), 
                             weight = model_1$wY) %>% 
    mutate(covariate = factor(covariate, levels = c("bipolar", "mdd", "schizo", 
                                                    "opioids", "anti_anxiety")))
  
  df_y_res_overall %>% 
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
  df_x_res_overall <- tibble(ensembl_gene_id = colnames(X.mat),
                             weight = model_1$wX) %>% 
    left_join(df_modules_filt) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = unique(.$ensembl_gene_id)))
  
  df_col <- df_modules_filt %>% 
    dplyr::select(module, color) %>% 
    distinct()
  colors <- df_col %>% pull(color)
  names(colors) <- df_col %>% pull(module)
  
  df_x_res_overall %>% 
    arrange(module) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = .$ensembl_gene_id)) %>% 
    
    ggplot(aes(x = ensembl_gene_id, y = weight, fill = module)) +
    geom_col() +
    scale_fill_manual(values = colors) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("X coefficient weights") +
    theme(legend.position = "bottom")
  
  

  
  
# Significant genes by bootstrapping ---------------------------------------------------

  
# PLOT Y RES
  
  df_y_res %>% 
    ggplot(aes(x = covariate, y = weight)) +
    geom_point(size = 5, aes(color = weight)) +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    scale_color_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
                          limits = c(-0.1, 0.1)) +
    ylim(c(-0.1, 0.1)) +
    xlab("") +
    ggtitle("Y coefficient weights") +
    theme(legend.position = "none")
  
  
# PLOT X RES
  df_x_res %>% 
    arrange(module) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = .$ensembl_gene_id)) %>% 
    
    ggplot(aes(x = ensembl_gene_id, y = weight, fill = module)) +
    geom_col() +
    scale_fill_manual(values = colors) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("X coefficient weights") +
    theme(legend.position = "none")
  

# PLOT MODULE ORDER
  df_x_res_mod <- df_x_res %>% 
    filter(weight != 0) %>% 
    count(module) %>% 
    left_join(df_modules_filt %>% 
                count(module) %>% 
                dplyr::rename("module_n" = "n")) %>% 
    mutate(perc_mod = n/module_n * 100) %>% 
    arrange(-perc_mod) %>% 
    mutate(module = factor(module, levels = .$module))
 
  df_x_res_mod %>% 
    ggplot(aes(x = perc_mod, y = module)) +
    geom_col(aes(fill = module, alpha = perc_mod)) +
    geom_text(aes(label = round(perc_mod, 1)), hjust = -0.05) +
    scale_fill_manual(values = colors) +
    # scale_fill_gradient(low = "white", high = "midnightblue",
    #                     limits = c(0, 25)) +
    labs(x = "percent of module in GRCCA significant gene list") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12))
  
# PLOT SIGNIFICANT GENES   
  df_x_res %>% 
    filter(module == 15 & weight != 0) %>% 
    arrange(-weight) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = .$gene_symbol)) %>% 

    ggplot(aes(x = weight, y = gene_symbol, 
               fill = abs(weight))) +
    geom_col() +
    scale_fill_gradient(low = "white", high = "midnightblue",
                         limits = c(0, 0.015)) +
    ylab("") +
    theme_classic() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12)
    ) +
    ggtitle("module 15 significant genes")
  
# CHECK DIFFERENTIAL EXPRESSION
  load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res
  
  grcca_genes <- df_x_res %>%
    filter(weight != 0) %>%
    pull(ensembl_gene_id) # significant
  
  df_akula_res %>% 
    filter(id %in% grcca_genes & padj < 0.05) %>% 
    left_join(df_modules_filt %>% dplyr::rename("id" = "ensembl_gene_id")) #%>% 
    write.csv("~/Downloads/akula_degs_grcca.csv")

    df_akula_res %>% 
      #filter(id %in% grcca_genes) %>% 
      filter(id %in% df_modules_filt$ensembl_gene_id & padj < 0.05) %>% 
      left_join(df_modules_filt %>% dplyr::rename("id" = "ensembl_gene_id")) %>% 
      count(module) %>% 
      arrange(-n) %>% 
      mutate(module = factor(module, levels = .$module)) %>% 
      
      ggplot(aes(x = n, y = module, fill = module)) +
      geom_col() +
      scale_fill_manual(values = colors) +
      theme(legend.position = "none")
    
  
# overlap with SCZ risk genes ---------------------------------------------


# data
  load("objects/df_hsapiens_genome.RDS") # df_hsapiens_genome
  df_scz_genome <- readxl::read_excel("data/eva/lago_bahn_2022_scz_genome.xlsx") %>% 
    row_to_names(1) %>% 
    clean_names
  
# hypergeometric test
  scz_genes <- df_scz_genome$risk_gene_hgnc %>% unique()
  
  gene_universe <- df_hsapiens_genome %>% 
    filter(gene_id %in% colnames(X.mat)) %>% 
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>% 
    pull(gene_name) %>% unique
  
  scz_genes_in_universe <- intersect(gene_universe, scz_genes) %>% unique
  
  grcca_genes <- df_x_res %>%
    filter(weight != 0) %>% 
    pull(gene_symbol) # 1219
  
  overlap_genes <- intersect(grcca_genes, scz_genes_in_universe)
  
  q <- length(overlap_genes) - 1
  m <- length(scz_genes_in_universe)
  n <- length(gene_universe) - m
  k <- length(grcca_genes)
  phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE) # p = 0.021
  
# plot overlap genes
  df_x_res %>% 
    filter(gene_symbol %in% overlap_genes) %>% 
    arrange(weight) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = unique(.$gene_symbol))) %>% 
    ggplot(aes(x = weight, y = gene_symbol, fill = module)) +
    geom_col() +
    scale_fill_manual(values = colors) +
    # scale_fill_gradient2(low = "midnightblue", mid = "white", high = "midnightblue",
    #                      limits = c(-0.02, 0.02)) +
    ylab("") +
    theme_classic() +
    theme(plot.title = element_text(size = 18),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 15),
          legend.position = "bottom"
    ) +
    ggtitle("SCZ druggable targets; GRCCA sig genes")
  
# chromosome region
  
  # get gene info from ensembl database
  hsapiens_ensembl <- useEnsembl(biomart = "genes",
                                 dataset = "hsapiens_gene_ensembl")

  df_chromosomes <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                          filters = "hgnc_symbol",
                          values = overlap_genes,
                          mart = hsapiens_ensembl) %>% 
    tibble() %>% 
    filter(!grepl("CHR", chromosome_name)) %>% 
    mutate(chromosome_name = as.numeric(chromosome_name)) %>% 
    arrange(chromosome_name, start_position) %>% 
    left_join(read_csv("data/chromosome_lengths.csv")) %>% 
    mutate(chromosome_name = factor(paste0("chr", chromosome_name), 
                                    levels = paste0("chr", unique(.$chromosome_name)))) 
  
  
  df_chromosomes %>% 
    left_join(df_modules_filt %>% 
                dplyr::rename("gene_id" = "ensembl_gene_id") %>% 
                left_join(df_hsapiens_genome %>% 
                            dplyr::select(gene_id, gene_name) %>% 
                            distinct() %>% 
                            filter(gene_name %in% overlap_genes)) %>% 
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
                  ymin = start_position, ymax = end_position, fill = module,
                  fill = hgnc_symbol)) +
    geom_label_repel(aes(y = start_position, label = hgnc_symbol, fill = module),
                     alpha = 0.8, min.segment.length = 0.1,
                     max.overlaps = 30) +
    scale_fill_manual(values = colors) +
    #facet_wrap(vars(chromosome_name), scales = "free") +
    coord_flip() +
    labs(x = "", y = "") +
    theme_void() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 12),
          legend.position = "none") +
    ggtitle("Genomic region of SCZ risk genes")
  
  
  


  
