
########################################################################################

# GENE LEVEL: Associate module PCs with covariates 

########################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/03.WGCNA/load_WGCNA_res.R"))


# setup ---------------------------------------------------------------

## LIBRARIES
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


## PLOT THEME
theme_set(theme_bw() +
            theme(plot.title = element_text(size = 12),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 10),
                  strip.text = element_text(size = 10),
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 10)))

# dx colors
dx_colors <- c("#0072B2", "#E69F00", "#009E73", "#9966FF")
names(dx_colors) <- c("Control", "BD", "MDD", "SCZ")

## SET DIRECTORIES
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "08Mar2024_GENES_qSVAgeSexRaceGC"

## IDENTIFY MODULE SET OF INTEREST
soft_power <- 3
minimum_size <- 40
tree_cut_height <- 0.98

# LOAD DATA OBJECTS 
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

# IDENTIFY MODULE SET OF INTEREST  
df_modules_filt <- df_modules %>% 
    filter(min_size == minimum_size & cut_height == tree_cut_height) %>% 
    unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module) %>% 
    
    # rename modules to specify gene-level
    mutate(module = paste0("geneM", module) %>% factor(levels = paste0("geneM", levels(module)))
    )

# Link dx number to group abbreviation
group <- c(0, 1, 2, 3)
names(group) <- c("BD", "Control", "MDD", "SCZ")

# DRUG MCA RESULTS
load(paste0(base_dir, "objects/drug_MCA_results.Rdata")) # df_ind_loadings, df_var_loadings

# DEFINE COVARIATE TYPES
technical_covariates <- c("source", "pmi", "ph", "pmi_confidence", "max_rin", "max_rine", "mapped_percent", "five_prime_three_prime_bias",
                          "rna_extraction_batch", "library_batch", "gc_percent")
biological_covariates <- c("age_death", "sex_at_birth", "race", "marital_status", "education", "manner_death", "suicide", "brain_weight",
                           "height", "weight", "bmi")
drug_covariates <- c("smoker", "nicotine_cotinine", "alcohol", "sedative_hypnotic_anxiolitics", "opioids", "cannabinoids",
                     "major_stimulants_cocaine_included", "minor_stimulants", "anticholinergics", "antidepressants",
                     "anti_epileptics", "anti_histamines", "antipsychotics", "mood_stabilizers", "non_psychiatric",
                     "other_psychotropic_drug", "benzos")

# ADD DRUG MCA DIMENSIONS TO COVARIATE DF
df_covariates_mca <- df_covariates_numeric %>% 
  left_join(
    df_ind_loadings %>% 
      dplyr::rename_at(vars(contains("dim")), ~ str_replace(.x, "dim", "MC")),
    by = join_by(sample)
  ) %>% 
  
  # remove drug covariates (since we have the MCs)
  dplyr::select(-any_of(drug_covariates))



# run PCA on each module --------------------------------------------------


# RUN PCA FOR EACH MODULE
df_mods_pca <- df_modules_filt %>% 
  
  # combin with sample expression information
  left_join(df_vsd_regress %>% 
              pivot_longer(2:ncol(.), names_to = "ensembl_gene_id", values_to = "resids"),
            by = join_by(ensembl_gene_id)
  ) %>% 
  
  # create transcript x sample expression matrix for each module
  pivot_wider(id_cols = c(ensembl_gene_id, module), 
              names_from = sample, 
              values_from = resids) %>% 
  group_by(module) %>% 
  nest() %>% 
  
  # run PCA on transcript x sample expression matrix
  mutate(pca = map(.x = data,
                   .f = ~ .x %>% 
                     as.data.frame %>%
                     column_to_rownames("ensembl_gene_id") %>% 
                     prcomp(center = TRUE, scale = TRUE)
  )) %>% 
  
  # combine PCA results
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
  
  # combin sample PC scores with sample covariate information
  left_join(df_covariates_mca %>% 
              pivot_longer(2:ncol(.), 
                           names_to = "covariate", 
                           values_to = "covariate_val"),
            by = "sample")

df_mods_pca


# quantify module-trait associations --------------------------------------


# CORRELATE EACH PC WITH ALL COVARIATES
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
  group_by(module) %>% 
  mutate(q_val = p.adjust(p_val, method = "fdr")) %>% 
  mutate(pc = str_remove(pc, "pc") %>% 
           as.numeric %>% 
           as.factor) %>% 
  arrange(p_val)

  # print significant associations
  df_mods_pc_covar_cor %>% 
    filter(q_val < 0.05) %>% 
    arrange(pc, module) %>% 
    print(n = nrow(.))

# WHAT COVARIATES RELATE TO SIGNIFICANT MODULE PCs?
df_mods_pc_covar_cor %>% 
  unnest(cols = c(data)) %>% 
  filter(proportion_of_variance > 0.02 & p_val < 0.05) %>% 
  ungroup %>% 
  dplyr::select(module, covariate) %>% 
  distinct() %>% 
  dplyr::count(covariate, sort = TRUE)

  # correlates with eigengene only
  df_mods_pc_covar_cor %>% 
    filter(pc == 1 & p_val < 0.05) %>% 
    ungroup %>% 
    dplyr::select(module, covariate) %>% 
    distinct() %>% 
    dplyr::count(covariate, sort = TRUE) %>% 
    print(n = nrow(.))
  
  df_mods_pc_covar_cor %>% 
    filter(pc == 1 & p_val < 0.05 & covariate %in% technical_covariates)


# PLOT

# set order
m_covariates <- df_covariates_mca %>% column_to_rownames("sample") %>% t()
covariate_order <- m_covariates[hclust(dist(m_covariates))$order,] %>% rownames

# plot
df_mods_pc_covar_cor %>% 
  unnest(cols = c(data)) %>% 
  mutate(covariate = factor(covariate, levels = covariate_order)) %>% 
  group_by(module) %>% 
  filter(proportion_of_variance > 0.02) %>% # filter(q_val < 0.05) %>% select(module, pc, covariate) %>% distinct
  mutate(pearsons_r = ifelse(p_val > 0.05, NA, pearsons_r)) %>% 
  
  ggplot(aes(x = pc, y = covariate)) +
  geom_tile(aes(fill = pearsons_r)) +
  scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC", na.value = "white",
                       limits = c(-1, 1)) +
  facet_wrap(~ module, scales = "free") +
  labs(x = "PC explained var > 2%", title = "Module PC-covariate correlations") +
  theme(#axis.text.y = element_blank(),
    #axis.text.x = element_text(size = 8),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10))




# Fig SXX ?? | extract eigenegene, plot example ----------------------------------

# GET EIGENGENES
df_kme <- df_mods_pca %>% 
  filter(pc == "pc1") %>% 
  pivot_wider(id_cols = c(module, sample, pc_score), 
              names_from = covariate, 
              values_from = covariate_val) %>% 
  dplyr::rename("kme" = "pc_score") %>% 
  left_join(enframe(group) %>% dplyr::rename("dx" = "value")) %>% 
  dplyr::select(-dx) %>% 
  dplyr::rename("dx" = "name") %>% 
  mutate(dx = factor(dx, levels = c("Control", "BD", "MDD", "SCZ"))) %>% 
  dplyr::select(module, sample, kme, dx, everything())

# EXAMPLE PLOT OF WHAT IS AN EIGENGENE
  df_mods_pca %>% 
    filter(module == 6) %>% 
    filter(pc %in% c("pc1", "pc2")) %>% 
    mutate(pc = paste0(toupper(pc), 
                       sep = "; proportion of variance explained = ", 
                       round(proportion_of_variance, 3))) %>% 
    dplyr::select(sample, pc, pc_score) %>% distinct %>% 
    pivot_wider(id_cols = sample, 
                names_from = pc,
                values_from = pc_score) %>% 
    
    mutate(dx = case_when(
               str_detect(sample, "control") ~ "Control",
               str_detect(sample, "bipolar") ~ "BD",
               str_detect(sample, "mdd") ~ "MDD",
               str_detect(sample, "schizo") ~ "SCZ"
             ) %>% 
             factor(levels = c("Control", "BD", "MDD", "SCZ"))) %>% 
    
    ggplot(aes(x = `PC1; proportion of variance explained = 0.151`,
               y = `PC2; proportion of variance explained = 0.065`,
               fill = dx)) +
    geom_point(shape = 21, size = 2.5, color = "black") +
    geom_vline(aes(xintercept = 0), color = "black", lty = 2) +
    geom_hline(aes(yintercept = 0), color = "black", lty = 2) +
    scale_fill_manual(values = dx_colors) +
    labs(title = bquote(bold("A")~" | Module 6 Principal Components Analysis")) +
    theme(legend.position = c(0.9, 0.15),
          legend.key.height = unit(0.4, "cm"))
  
  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0("~/Documents/PhD/manuscripts/BiolPsych2024/figures/supplement/FigS3A", .x),
                  width = 6, height = 5)
  )
  
  
  
# S6B | relate module eigengene to diagnosis ------------------------------

## EXTRACT PC1 AS EIGENGENE
  df_kme <- df_mods_pca %>% 
      
      # pull eigengene value for each sample
      filter(pc == "pc1") %>% 
      pivot_wider(id_cols = c(module, sample, pc_score), 
                  names_from = covariate, 
                  values_from = covariate_val) %>% 
      dplyr::rename("kme" = "pc_score") %>% 
      mutate(dx = case_when(
          
          dx == 0 ~ "BD",
          dx == 1 ~ "Control",
          dx == 2 ~ "MDD",
          dx == 3 ~ "SCZ"
          
      ) %>% 
          factor(levels = c("Control", "BD", "MDD", "SCZ")),
      .before = 4)
  
  
# LINEAR REGRESSION OF MES AND DIAGNOSIS
  df_lm_dx <- df_kme %>% 
      group_by(module) %>% 
      nest() %>% 
      mutate(
          data = map(
              .x = data,
              .f = ~ .x %>% 
                  mutate(kme_resids = lm(kme ~ brain_weight) %>% residuals + median(kme))
          ),
          lm = map(
              .x = data,
              .f = ~ lm(kme_resids ~ dx + MC1 + MC2 + MC3 + MC4 + MC5 + MC6 + MC7 + MC8, data = .x)
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
  
## CALCULATE RESIDUALS FOR PLOTTING
  df_lm_dx_plot <- df_lm_dx %>% 
      mutate(
          data = map(
              .x = data,
              .f = ~ .x %>% 
                  mutate(kme_resids = lm(kme ~ brain_weight + MC1 + MC2 + MC3 + MC4 + MC5 + MC6 + MC7 + MC8) %>% residuals + median(kme))
          )
      ) %>% 
      filter(p_value < 0.05 & term %in% names(dx_colors)) %>% 
      unnest(cols = c(data)) %>% 
      dplyr::select(module, term, dx, kme_resids, p_value, p_adj) %>% 
      
      # add star for significance
      mutate(significant = ifelse(term == dx & p_value < 0.05, "*", ""))
  

# PLOT SIGNIFICANT MODULE EIGENGENE-DX ASSOCIATIONS
  df_lm_dx_plot %>% 
      
      # plot
      ggplot(aes(x = kme_resids, y = dx, color = dx)) +
      geom_point(aes(fill = dx), color = "black", shape = 21, alpha = 0.5,
                 position = position_jitter(height = 0.25)) +
      geom_boxplot(fill = "transparent", outlier.shape = NA) +
      geom_text(aes(label = significant), x = -0.15, color = "black",
                size = 7, vjust = 0.75) +
      scale_color_manual(values = dx_colors) +
      scale_fill_manual(values = dx_colors) +
      facet_wrap(vars(module)) +
      labs(y = "", x = "Module eigengene residuals",
           caption = "kME residuals ~ diagnosis + MC1 + MC2 + MC3 + MC4 + MC5 + MC6 + MC7 + MC8",
           title = "Modules with significant case-control module eigengene differences") +
      theme(legend.position = "none")
  
  # save to project dir
  map(
      .x = c(".png", ".pdf"),
      .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/genes/S6A.module_kME_diagnosis", .x),
                    width = 8, height = 4)
  )
 
# PLOT OTHER SIGNIFICANT ASSOCIATIONS
  df_lm_dx %>% 
      filter(p_value < 0.05) %>% 
      unnest(cols = c(data)) %>% 
      
      filter(module == 9 & term == "MC7") %>% 
      ggplot(aes(x = MC7, y = kme_resids)) +
      geom_point() +
      geom_smooth(method = "lm") +
      stat_cor()
  
   
# EXPORT TABLE
  df_lm_dx %>% 
      filter(p_value < 0.05) %>% 
      dplyr::select(module, term, estimate, std_error, p_value, p_adj) %>% 
      bind_rows(tibble(module = "kME (corrected for brain weight) ~ diagnosis + MC1 + MC2 + MC3 + MC4 + MC5 + MC6 + MC7 + MC8"), .) %>% 
      write_xlsx(paste0(base_dir, "outputs/tables/for_manuscript/TableS2_GENES_module_eigenegene_res.xlsx"))
  

# Fig SXX | Module PC - covariate correlations --------------------------------------
  
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
    mutate(p_adj = p.adjust(p_val, method = "fdr")) %>% 
    mutate(pc = str_remove(pc, "pc") %>% as.numeric %>% as.factor)
  
  df_mods_pc_covar_cor %>% 
    unnest(cols = c(data)) %>% 
    filter(proportion_of_variance > 0.02 & p_val < 0.05) %>% 
    dplyr::select(module, pc, covariate, pearsons_r, p_val, p_adj) %>% 
    distinct() %>% 
    arrange(module, pc) %>% 
    print(n = nrow(.))
  
# PLOT
  technical_covars <- c("rna_extraction_batch", "pmi_confidence", "max_rine", "max_rin", "gc_percent",
                        "source", "pmi", "ph", "mapped_percent", "library_batch", "five_prime_three_prime_bias")
  
  df_mods_pc_covar_cor %>% filter(module == 6) %>% 
    unnest(cols = c(data)) %>% 
    group_by(module) %>% 
    filter(proportion_of_variance > 0.01) %>% # filter(q_val < 0.05) %>% select(module, pc, covariate) %>% distinct
    mutate(pearsons_r = ifelse(p_val > 0.05, NA, pearsons_r)) %>% 
    mutate(labels = case_when(
      p_val < 0.05 & p_adj > 0.05 ~ "*",
      p_adj < 0.05 ~ "**"
    )) %>% 
    filter(!(covariate %in%  technical_covars) &
             !str_detect(covariate, "qSV")) %>% 
    mutate(covariate = str_replace_all(covariate, "_", " ")) %>% 
    dplyr::rename("r" = "pearsons_r") %>% 
    
    # plot
    ggplot(aes(x = pc, y = covariate)) +
    geom_tile(aes(fill = r)) +
    geom_text(aes(label = labels), color = "black", size = 5, vjust = 0.75,
              check_overlap = TRUE) +
    scale_fill_gradient2(high = "#B2182B", mid = "white", low = "#2166AC", na.value = "white",
                         limits = c(-1, 1)) +
    facet_wrap(~ module, scales = "free_x", nrow = 3) +
    labs(x = "PC explained variance > 1%", y = "",
         title = bquote(bold("C")~" | Module PC - covariate correlations"),
         caption = "Many PCs beyond PC1 contribute to variation in data") +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          strip.text = element_text(size = 10))

  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0("~/Documents/PhD/manuscripts/BiolPsych2024/figures/supplement/FigS3C", .x),
                  width = 8, height = 8)
  )
  

# Fig SXX | Significant role of drugs in kME (leave for now) ----------------------------------------------------------------------

  drug_covariates <- c("smoker", "nicotine_cotinine", "alcohol", "sedative_hypnotic_anxiolitics", "opioids", "cannabinoids",
                       "major_stimulants_cocaine_included", "minor_stimulants", "anticholinergics", "antidepressants",
                       "anti_epileptics", "anti_histamines", "antipsychotics", "mood_stabilizers", "non_psychiatric",
                       "other_psychotropic_drug", "benzos")
  
  
# LINEAR REGRESSION OF MES AND DRUGS
  df_lm_drugs <- df_kme %>% 
    dplyr::select(module, kme, all_of(drug_covariates)) %>% 
    group_by(module) %>% 
    nest() %>% 
    mutate(
      
      lm = map(
        .x = data,
        .f = ~ lm(kme ~ ., data = .x)
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
    mutate(FDR = p.adjust(p_value, method = "BH"),
           significance = case_when(
             FDR < 0.05 ~ "**",
             FDR > 0.05 & p_value < 0.05 ~ "*",
             p_value > 0.05 ~ ""
           )
    )
  
# PLOT KME DISTRIBUTION ACROSS DRUGS
  df_kme %>%
    left_join(df_lm_drugs, by = join_by(module, dx)) %>% 
    dplyr::select(module, kme, dx, significance) %>% 
    distinct() %>% #filter(significance != "") %>% 
    
    ggplot(aes(x = kme, y = dx, color = dx)) +
    geom_violin() +
    stat_summary(fun = "median", geom = "crossbar", width = 1) +
    geom_point(position = position_jitter(height = 0.25), shape = 1, alpha = 0.7) +
    geom_text(aes(label = significance), x = 0.17, vjust = 0.75, 
              size = 8, color = "black", check_overlap = TRUE) +
    geom_vline(aes(xintercept = 0), color = "black", lty = 2) +
    scale_color_manual(values = dx_colors) +
    facet_wrap(vars(module)) +
    labs(y = "", 
         title = bquote(bold("B")~" | Module eigenegene distributions across diagnostic groups")) +
    theme(legend.position = "none")
  
  map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0("~/Documents/PhD/manuscripts/BiolPsych2024/figures/supplement/FigS3B", .x),
                  width = 10, height = 8)
  )
  