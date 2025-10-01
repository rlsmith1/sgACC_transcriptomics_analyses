
# ---------------------------------------------------------------------------- #
# Rerun toxicology MCA treating missing variables by filling in the
# individual-level mode (instead of marking as absent)
# ---------------------------------------------------------------------------- #


# Setup -------------------------------------------------------------------


### SET DIRECTORIES ###

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/response_to_reviewers/round2/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")


### LOAD LIBRARIES ###

library(tidyverse)
library(DESeq2)
library(MASS)

library(ggh4x)
library(pals)


### LOAD DATA ###

## Cleaned covariates
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in 00.clean_covariates.R)


## Define toxicological covariates
drug_covariates <- c(
    "smoker", "nicotine_cotinine", "alcohol", "opioids", "cannabinoids", "beta_blockers",
    "major_stimulants_cocaine_included", "minor_stimulants", "anticholinergics", "antidepressants",
    "anti_epileptics", "anti_histamines", "antipsychotics", "mood_stabilizers", "sedative_hypnotic_anxiolitics", "benzos",
    "other_psychotropic_drug", "non_psychiatric"
)
technical_covariates <- c(
    "mapped_percent", "gc_percent", "five_prime_three_prime_bias", "rin_acsg", "rna_extraction_batch",
    "library_batch", "pmi", "pmi_confidence", "source", "ph", "max_rin", "max_rine"
)
biological_covariates <- c(
    "age_death", "sex_at_birth", "race", "bmi", "height", "weight",
    "marital_status", "manner_death", "education", "suicide", "brain_weight"
)


## Define diagnosis colors (for plotting)
## Dx colors
dx_colors <- c(
    "Control" = "#0072B2", 
    "BD" = "#E69F00", 
    "MDD" = "#009E73", 
    "SCZ" = "#9966FF"
)



# Filter data and address missing values ----------------------------------------------------------


## Select drugs for analysis
df_drugs <- df_covariates_numeric %>% 
    ungroup %>% 
    dplyr::select(sample, all_of(drug_covariates)) %>% 
    distinct


## For each variable, calculate the percent that is missing values
exclude_drugs <- df_drugs %>%
    dplyr::select(-sample) %>% 
    summarise(across(everything(), ~ (sum(is.na(.)) / n()) * 100)) %>% 
    pivot_longer(1:ncol(.), names_to = "drug", values_to = "perc_missing") %>% 
    
    # include only drug covariates that retain at least 66.6% of the original observations
    dplyr::filter(perc_missing > 33.3) %>% 
    pull(drug) # this remove beta blockers, anti-histamines, and mood stabilizers
include_drugs <- df_drugs %>% colnames %>% setdiff(c("sample", exclude_drugs))


## Write helpfer function to replace NAs with the mode
get_mode <- function(x) {
    # Remove missing values first
    x_clean <- x[!is.na(x)]
    
    if (length(x_clean) == 0) {
        # If the row is all NA, return NA (cannot impute)
        return(NA)
    }
    
    # Calculate mode
    ux <- unique(x_clean)
    mode_val <- ux[which.max(tabulate(match(x_clean, ux)))]
    
    return(mode_val)
}


## For all missing values, fill in with the individual-level mode
## Rationale: clinically, if an individual is taking one drug they are more likely to be taking another
df_mca_drugs <- df_drugs %>%
    
    # remove drugs with too much (>33.3% missing data)
    dplyr::select(all_of(include_drugs)) %>% 
    
    # Create a new column with the row-wise mode for each subject
    rowwise() %>% 
    mutate(row_mode = get_mode(c_across(all_of(include_drugs)))) %>%
    ungroup %>%
    
    # Impute NAs in the drug columns using
    mutate(across(all_of(include_drugs), ~ ifelse(is.na(.), row_mode, .))) %>% 
    
    # Remove the temporary row mode column
    dplyr::select(-row_mode) %>% 

    # Set all drug column variables to factors
    mutate_all( ~ factor(.x)) %>% 
    
    # Add sample information back
    mutate(sample = df_drugs$sample, .before = 1) %>% 
    
    # Convert to dataframe object for MCA algorithm
    column_to_rownames("sample")



# Run MCA -----------------------------------------------------------------

## Define number of factors to calculate
n_factors <- length(include_drugs) # 15 drugs to include


## Run MCA
set.seed(20250929)
mca.res <- mca(df_mca_drugs,
               nf = n_factors
)

## Calculate eigenvalue of each dimension
df_eig <- tibble(
    dim = 1:n_factors,
    var = (mca.res$d^2)*100,
    cum_var = cumsum(var)
)

## Extract variable loadings
df_var_loadings <- mca.res$cs %>% 
    as.data.frame %>% 
    dplyr::rename_all( ~ paste0("dim", .x)) %>% 
    rownames_to_column("variable") %>% 
    as_tibble() %>% 
    separate(variable, into = c("drug", "use"), sep = "\\.") %>% 
    mutate(use = ifelse(use == 1, "yes", "no"))

## Extract individual loadings on each dimension
df_ind_loadings <- mca.res$rs %>%
    as.data.frame %>%
    dplyr::rename_all( ~ paste0("dim", .x)) %>%
    rownames_to_column("sample") %>%
    as_tibble()



# Correlate each covariate with each MCA dimension -----------------------------------------


## Dummy code diagnosis variable for correlation
all(df_covariates_numeric$age_death == df_covariates$age_death) # check row order
df_covariates_dxDummy <- df_covariates_numeric %>% 
    mutate(dx = df_covariates$dx) %>% 
    mutate(
        Control = ifelse(dx == "Control", 1, 0),
        BD = ifelse(dx == "BD", 1, 0),
        MDD = ifelse(dx == "MDD", 1, 0),
        SCZ = ifelse(dx == "SCZ", 1, 0)
    ) %>% 
    dplyr::select(-dx)


## Correlate each covariate with each MCA dimension
df_mca_covar_cor <- expand_grid(
    covariate = df_covariates_dxDummy %>% dplyr::select(-sample) %>% colnames,
    dimension = df_ind_loadings %>% dplyr::select(-sample) %>% colnames
) %>% 
    mutate(
        pearsons_r = map2(
            .x = covariate,
            .y = dimension,
            .f = ~ cor.test(df_covariates_dxDummy[[.x]], df_ind_loadings[[.y]])$estimate
        ),
        p_val = map2(
            .x = covariate,
            .y = dimension,
            .f = ~ cor.test(df_covariates_dxDummy[[.x]], df_ind_loadings[[.y]])$p.value
        )
    ) %>% 
    unnest(cols = c(pearsons_r, p_val)) %>% 
    mutate(p_adj = p.adjust(p_val, method = "fdr"))



# Save for sensitivity analysis figures -------------------------------------------


## Save for figures
save(df_eig, df_var_loadings, df_ind_loadings, df_mca_covar_cor,
     file = paste0(analysis_objects_dir, "R2-drug_MCA_results_modeImpute.Rdata")
)



# FIG S3A v2 -----------------------------------------------------------------


## Create numeric matrix of imputed toxicology data
df_drugs_impute <- df_mca_drugs %>%
    
    # Set rownames as sample ID column
    rownames_to_column("sample") %>% 
    as_tibble() %>% 
    
    # Add diagnosis label to group by; will cluster within in each diagnosis group
    mutate(dx = str_remove(sample, ".*_"), .before = 2) %>% 
    
    # Make sure all columns are numeric (not factor)
    mutate_at(vars(all_of(include_drugs)), ~ as.numeric(.x) - 1) 



## Cluster samples and drugs using hclust to set order for plot
df_drugs_cluster <- df_drugs_impute %>% 
    
    # Run clustering algorithm within each diagnosis group
    group_by(dx) %>% 
    nest() %>%
    
    # Run hclust within each diagnosis group to cluster
    mutate(
        sample_hclust = map(
            .x = data,
            .f = ~ hclust(dist(.x %>% column_to_rownames("sample")))
        ),
        drugs_hclust = map(
            .x = data,
            .f = ~ hclust(dist(.x %>% column_to_rownames("sample") %>% t()))
        )
    )


## Identify order for samples based on hclust algorithm
samples_order <- df_drugs_cluster %>%
    
    # Pull sample order 
    mutate(order = map(.x = sample_hclust, .f = ~ .x$order)) %>% 
    unnest(cols = c(data, order)) %>% 
    dplyr::select(dx, order, sample) %>%
    
    # Set diagnosis as numeric variable to order
    mutate(dx = case_when(
        str_detect(sample, "control") ~ 1,
        str_detect(sample, "bipolar") ~ 2,
        str_detect(sample, "mdd") ~ 3,
        str_detect(sample, "schizo") ~ 4
    )) %>%
    
    # Arrange by diagnosis and then within-diagnosis hclust order
    arrange(dx, order) %>% 
    pull(sample)


## Identify order for drugs based on hclust algorithm
drugs_order <- hclust(dist(df_drugs_impute %>% dplyr::select(-dx) %>% column_to_rownames("sample") %>% t()))$order
drugs_order <- tibble(
    drug = df_drugs_impute %>% dplyr::select(-sample, -dx) %>% colnames,
    order = drugs_order
) %>% 
    arrange(order) %>% 
    pull(drug)


## Format tibble with arrange rows & columns
df_figSXXa <- df_drugs_impute %>% 
    
    # Pivot longer for tile plot
    pivot_longer(3:ncol(.), names_to = "drug", values_to = "value") %>% 
    
    # Add labels to indicate what 1 and 0 mean (unknown should not be present here since we imputed)
    mutate(
        
        value = case_when(
            value == 1 ~ "detected",
            value == 0 ~ "not detected",
            is.na(value) ~ "unknown"
        ),
        
        # Set sample order based on clustering
        sample = factor(sample, levels = samples_order),
        
        # Set drugs order based on clustering
        drug = factor(drug, levels = drugs_order),
        
        # Relabel diagnostic groups and set factor levels
        dx = case_when(
            dx == "bipolar" ~ "BD",
            dx == "control" ~ "Control",
            dx == "mdd" ~ "MDD",
            dx == "schizo" ~ "SCZ"
        ) %>% 
            factor(levels = names(dx_colors)),
        
    ) %>%

    # Add fill colors to indicate diagnostic group
    left_join(
        enframe(dx_colors, name = "dx", value = "dx_fill") %>%
            mutate(dx = factor(dx, levels = names(dx_colors)))
    )


## Set matrix colors for plot
matrix_colors <- c(
    "detected" = "white",
    "not detected" = "black",
    "unknown" = "gray"
)


## Count sample size within each diagnotic group (for plotting)
dx_sample_size <- df_figSXXa %>% dplyr::select(sample, dx) %>% distinct %>% count(dx) %>% deframe


## Plot!
df_figSXXa %>%
    
    # Plot layers
    ggplot(aes(x = sample)) +
    geom_tile(aes(y = drug, fill = value)) + # drug tiles
    geom_tile(aes(y = 0, fill = I(dx_fill))) + # sample tiles

    facet_wrap(vars(dx), nrow = 1, scales = "free_x", strip.position = "bottom") +
    force_panelsizes(cols = dx_sample_size) +
    
    # NEW: add diagnostic labels to samples (instead of legend)
    geom_text(aes(label = dx, x = 25, y = 0), color = "white", size = 3.0, check_overlap = TRUE) +

    # Aesthetics
    scale_fill_manual(values = matrix_colors, na.value = "gray") +
    guides(fill = guide_legend(title = NULL)) +
    labs(y = NULL, title = "A | Imputed sample toxicology data") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = c(-0.45, 0.90),
          legend.key.size = unit(0.2, "cm"),
          legend.text = element_text(size = 8),
          legend.box.background = element_rect(color = "black"),
          legend.box.margin = margin(t = 5, b = 5, l = 5, r = 5),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          axis.line = element_blank()
    )

## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C1A.sample_toxicology_matrix", .x), width = 6.5, height = 2.25)
)



# FIG S3B v2 -----------------------------------------------------------------

## MCA eigenvalues
df_eig %>% 
    mutate(dim = factor(dim, levels = rev(.$dim)),
           label = ifelse(cum_var < 80, paste0(round(cum_var, 1), "%"), NA_character_)
    ) %>% 
    
    # Plot
    ggplot(aes(y = dim)) +
    geom_col(aes(x = var, fill = var), color = "black") +
    geom_text(aes(x = var, label = label),
              size = 2.5, hjust = -0.1) +
    geom_vline(xintercept = 5, color = "maroon") +
    
    # Plot aesthetics
    scale_fill_gradient(low = "white", high = "black", guide = "none") +
    labs(x = "% variance explained", y = "Dimension",
         title = "B | Eigenvalues"
    )

## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C1B.MCA_eigenvalues", .x), width = 2.0, height = 2.5)
)



# FIG S3C v2 -----------------------------------------------------------------


## MCA variable loadings

## Prep data
figSXXC <- df_var_loadings %>% 
    pivot_longer(3:ncol(.), names_to = "dimension", values_to = "score") %>% 
    filter(use == "yes") %>% 
    mutate(
        drug = factor(drug, levels = rev(drugs_order)),
        dimension = str_replace(dimension, "dim", "") %>% factor(levels = 15:1),
        significant = ifelse(abs(score) > 0.015, "yes", "no")
    ) %>% 
    arrange(abs(score))

## Set color scale min and max
scale_max <- figSXXC %>% pull(score) %>% abs %>% max


## Plot 
figSXXC %>%
    
    # layout
    ggplot(aes(x = drug, y = dimension)) +
    geom_point(aes(fill = score, color = significant, size = abs(score)), shape = 21) +

    # aesthetics
    scale_fill_gradientn(colors = gene_weight_color_scale, limits = c(-scale_max, scale_max), name = "Loading") +
    scale_color_manual(values = c("yes" = "black", "no" = "transparent"), guide = "none") +
    scale_size_continuous(limits = c(0, scale_max), range = c(0.1, 5), guide = "none") +
    guides(fill = guide_colorbar(title = "Loading", title.position = "top", title.hjust = 0.5)) +
    
    coord_cartesian(clip = "off") +
    labs(y = NULL, x = NULL, 
         title = "C | Variable loadings") +
    theme(
        legend.margin = margin(l = -15),  
        legend.key.width = unit(0.15, "cm"),
        legend.key.height = unit(0.60, "cm"),
        legend.title = element_text(hjust = 0.5, size = 7),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 8.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
    )

## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C1C.MCA_variable_loadings", .x), width = 4.5, height = 3.3)
)



# FIG S3D v2 -----------------------------------------------------------------


load(paste0(analysis_objects_dir, "R2-drug_MCA_results_modeImpute.Rdata"))


## MCA covariate-dimension correlations

## Prep
df_figSXXD <- df_mca_covar_cor %>% 
    
    # Remove technical covariates & dx, as well as drugs that did not pass NA filtering
    # Only interested in drugs and biological/demographic covariates
    filter(covariate %in% c(names(dx_colors), biological_covariates, include_drugs, technical_covariates)) %>% 
    
    # Format for plotting
    mutate(
        
        dimension = str_remove(dimension, "dim") %>% as.numeric %>% factor(levels = 15:1),
        significant = ifelse(p_adj < 0.05, "yes", "no"),
        type = case_when(
            covariate %in% technical_covariates ~ "Technical",
            covariate %in% biological_covariates ~ "Demographic",
            covariate %in% drug_covariates ~ "Toxicological",
            TRUE ~ "Diagnosis"
        ) %>% 
            factor(levels = c("Diagnosis", "Demographic", "Technical", "Toxicological")),
        covariate = factor(covariate, levels = c(names(dx_colors), biological_covariates, technical_covariates, rev(drugs_order)))
        
    ) %>% 
    arrange(abs(pearsons_r))
    
## Plot
df_figSXXD %>% 
    
    # Layout
    ggplot(aes(x = covariate, y = dimension)) +
    geom_point(aes(fill = pearsons_r, size = abs(pearsons_r), color = significant), shape = 21) +

    facet_wrap(vars(type), scales = "free_x", nrow = 1) +
    force_panelsizes(cols = c(length(biological_covariates), length(drug_covariates))) +
    
    # Aesthetics
    scale_fill_gradientn(colors = rev(brewer.rdbu(100)), limits = c(-1.0, 1.0), name = "r") +
    scale_color_manual(values = c("yes" = "black", "no" = "transparent"), guide = "none") + #name = "FDR<0.05") +
    scale_size_continuous(limits = c(0, 1.0), range = c(0.1, 5), guide = "none") +
    
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Dimension", 
         title = "D | Dimension correlations with known covariates") +
    theme(
        legend.margin = margin(l = -15),  
        legend.key.width = unit(0.15, "cm"),
        legend.key.height = unit(0.60, "cm"),
        legend.title = element_text(hjust = 0.5, size = 7),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 8.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
    )


## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C1D.MCA_dimension_covariate_cor", .x), width = 6.5, height = 3.4)
)




# Fig SXX: Covariate correlation plot including diagnosis and technical variables -----------------------------------------------------------------


## Plot
df_figSXXD %>% 
    
    # Layout
    ggplot(aes(x = covariate, y = dimension)) +
    geom_point(aes(fill = pearsons_r, size = abs(pearsons_r), color = significant), shape = 21) +
    
    facet_wrap(vars(type), scales = "free_x", nrow = 1) +
    force_panelsizes(cols = c(
        length(dx_colors), 
        length(biological_covariates),
        length(technical_covariates),
        length(drug_covariates)
    )
    ) +
    
    # Aesthetics
    scale_fill_gradientn(colors = rev(brewer.rdbu(100)), limits = c(-1.0, 1.0), name = "r") +
    scale_color_manual(values = c("yes" = "black", "no" = "transparent"), guide = "none") + #name = "FDR<0.05") +
    scale_size_continuous(limits = c(0, 1.0), range = c(0.1, 4), guide = "none") +
    
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Dimension", 
         title = "MCA dimension correlations with known covariates") +
    theme(
        legend.margin = margin(l = -15),  
        legend.key.width = unit(0.15, "cm"),
        legend.key.height = unit(0.60, "cm"),
        legend.title = element_text(hjust = 0.5, size = 7),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 8),
        axis.text.y = element_text(size = 8.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
    )


## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C1.MCA_dimension_covariate_cor_ALL", .x), width = 6.5, height = 3.4)
)


