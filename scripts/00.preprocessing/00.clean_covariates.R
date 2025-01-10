
######################################################################################

# Tidy covariate data for downstream analyses

######################################################################################


# libraries ---------------------------------------------------------------

library(tidyverse)
library(readxl)


# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"

# complete covariates from Ajeet
df_complete_covariates <- read_xlsx(paste0(base_dir, "data/14Nov2023_all_covariates.xlsx"))

# other covariates from Nirmala (includes additional technical covariates)
df_nirmala_covariates <- read.table(paste0(base_dir, "data/185samples_multinomialDx_46covariates_10evecs.txt"), 
                                    sep = "\t", header = TRUE) %>% 
  dplyr::select(-contains("evec")) %>% 
  as_tibble() %>% 
  mutate_if(is.character, ~ ifelse(.x == "", NA_character_, .x)) %>% 
  clean_names %>% 
  dplyr::rename("rna_extraction_batch" = "rn_aextraction_batch", 
                "gc_percent" = "g_cpercent",
                "dx" = "multi_nomial_dx") %>% 
  dplyr::select(-contains("_1")) %>% # remove duplicate columns
  mutate(sample = tolower(sample))

# define covariate types
technical_covariates <- c("source", "pmi", "ph", "pmi_confidence", "max_rin", "max_rine", "mapped_percent", "five_prime_three_prime_bias",
                          "rna_extraction_batch", "library_batch", "gc_percent")
biological_covariates <- c("age_death", "sex_at_birth", "race", "marital_status", "education", "manner_death", "suicide", "brain_weight",
                           "height", "weight", "bmi")
drug_covariates <- c("smoker", "nicotine_cotinine", "alcohol", "sedative_hypnotic_anxiolitics", "opioids", "cannabinoids",
                     "major_stimulants_cocaine_included", "minor_stimulants", "anticholinergics", "antidepressants",
                     "anti_epileptics", "anti_histamines", "antipsychotics", "mood_stabilizers", "non_psychiatric",
                     "other_psychotropic_drug", "benzos")


# format data -------------------------------------------------------------

df_complete_covariates2 <- df_complete_covariates %>% 
  clean_names %>%
  dplyr::select(-c(brain_number_5, 
                   verified, 
                   inclusion, 
                   avaialble_short_read_rna_seq, 
                   pseudo_guid, 
                   primary_dx,
                   gender_at_death,
                   sexual_orientation,
                   ethnicity,
                   contains("axi"),
                   cause_death,
                   final_neuropath_diagnosis,
                   inhalants_zinhal, # (only 1 non-NA value)
                   hallucinogens_zhall # (everyone is negative)
  )
  ) %>% 
  dplyr::rename("brain_id" = "brain_number_1", 
                "ph" = "p_h",
                "max_rine" = "max_ri_ne") %>% 
  dplyr::rename_all(~str_remove(.x, "_z.*")) %>% 
  mutate(sample = paste0(brain_id, sep = "_", tolower(dx)), .before = 1) %>% 
  mutate(suicide = ifelse(manner_death == "Suicide", 1, 0), .before = brain_weight) %>% 
  mutate(dx = case_when(
    dx == "Schizo" ~ "SCZ",
    dx == "Bipolar" ~ "BD",
    TRUE ~ dx
  ))

shared_variables <- intersect(colnames(df_complete_covariates2), colnames(df_nirmala_covariates))
colnames(df_nirmala_covariates)[!(colnames(df_nirmala_covariates) %in% shared_variables)]

# add technical covariates from Nirmala's df
df_covariates <- df_complete_covariates2 %>% 
  left_join(
    df_nirmala_covariates %>% 
      dplyr::select(sample, 
                    mapped_percent, 
                    five_prime_three_prime_bias, 
                    rna_extraction_batch, 
                    library_batch, 
                    gc_percent,
                    benzos),
    by = join_by(sample)
  )

# convert each covariate to numeric
df_covariates_numeric <- df_covariates %>% 
  dplyr::select(-brain_id) %>% 
  mutate_if(is.character, ~ ifelse(.x == "NULL", NA_character_, .x)) %>% 
  mutate_if(is.character, ~ as.factor(.x) %>% as.numeric - 1) %>% 
  mutate(sample = df_covariates$sample)

# SAVE OBJECTS
save(df_covariates, df_covariates_numeric, file = paste0(base_dir, "objects/14Nov2023_covariates.Rdata"))

## TABLE OF N
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata"))
df_covariates %>% count(sex_at_birth)
df_covariates %>% count(race)


# PLOT: correlate covariates with each other ------------------------------------

load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata"))

# calculate correlations
df_covariate_correlations <- expand_grid(
  cov1 = df_covariates_numeric %>% dplyr::select(-sample) %>% colnames,
  cov2 = df_covariates_numeric %>% dplyr::select(-sample) %>% colnames
) %>% 
  mutate(
    pearsons_r = map2(
      .x = cov1,
      .y = cov2,
      .f = ~ cor.test(df_covariates_numeric[[.x]], df_covariates_numeric[[.y]])$estimate
    ),
    p_value = map2(
      .x = cov1,
      .y = cov2,
      .f = ~ cor.test(df_covariates_numeric[[.x]], df_covariates_numeric[[.y]])$p.value
    )
  ) %>% 
  unnest(cols = c(pearsons_r, p_value)) %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

# set order
m_covariates <- df_covariates_numeric %>% column_to_rownames("sample") %>% t()
covariate_order <- m_covariates[hclust(dist(m_covariates))$order,] %>% rownames

# PLOT
df_covariate_correlations %>% 
  mutate(cov1 = factor(cov1, levels = covariate_order),
         cov2 = factor(cov2, levels = covariate_order)) %>% 
  mutate(pearsons_r = ifelse(cov1 == cov2, NA_integer_, pearsons_r),
         significant = ifelse(p_adj < 0.05 & cov1 != cov2, "yes", "no")) %>% 
  
  ggplot(aes(x = cov1, y = cov2)) +
  geom_tile(aes(fill = pearsons_r, color = significant), linewidth = 0.25) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       na.value = "white") +
  scale_color_manual(values = c("yes" = "black", "no" = "transparent"), guide = "none") +
  guides(fill = guide_colorbar(title = "Pearson's r")) +
  labs(x = "", y = "", title = "Covariate correlations",
       caption = "Box outline indicates FDR < 0.05") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.key.height = unit(2, "cm"))

# SAVE
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/covariates/covariate_correlations", .x),
                width = 10, height = 8)
)

