
########################################################################

# GENE LEVEL: Describe functional enrichment of modules implicated
# by univariate analysis and how they overlap with published DEG results

########################################################################


# Setup -------------------------------------------------------------------

## Set directories
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/WGCNA_res/")

## Load data and functions
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/load_WGCNA_res.R"))
source(paste0(base_dir, "scripts/full_analysis_scripts/functions/module_hypergeometric_overlap.R")) # f_module_hypergeometric

## Link dx to group number (in the numeric covariates tibble)
group <- c(0, 1, 2, 3)
names(group) <- c("BD", "Control", "MDD", "SCZ")


# Load data ---------------------------------------------------------------

### Study data ###

## module eigengene analysis results
load(paste0(base_dir, "objects/kME_analysis_res.RDS")) # df_lm_dx
sig_modules <- df_lm_dx %>% 
    dplyr::filter(p_value < 0.05 & term == "SCZ") %>% 
    pull(module)

## current analysis DEG results
load(paste0(base_dir, "objects/DE_results.RDS")) # df_de_res

### Benchmarking ###

## Akula DEGs
load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res

## Consensus DE data
df_consensus_degs <- readxl::read_excel(paste0(base_dir, "data/merikangas_2022_consensus_DEGs.xlsx")) %>% 
    left_join(df_ensembl_to_symbol)

## GWAS DATA
df_scz_gwas <- readxl::read_excel(paste0(base_dir, "data/trubetskoy_2023_scz_gwas.xlsx"),
                                  sheet = 3
)

# sheet 2 = 95% credible set
# sheet 3 = 95% credible set k <= 3.5
# sheet 4 = 95% credible set <= 5 SNPs
# sheet 5 = 95% credible set <= 5 PIP

## Rare variants
df_exome_urvs <- readxl::read_excel(paste0(base_dir, "data/singh_2022_exome_URVs.xlsx")) %>% 
    left_join(df_ensembl_to_symbol)



# Current study DEGs ------------------------------------------------------

### HYPERGEOMETRIC (overall enrichment) ###

# Pull DEG IDs
my_degs <- df_de_res %>% dplyr::filter(pvalue < 0.05) %>% pull(ensembl_gene_id)
length(my_degs) # n = 1747
my_degs %in% ensembl_ids %>% sum # n = 1747 in our universe

# Run hypergeometric test for each module
l_my_degs_hypergeometric <- map(
    .x = sig_modules,
    .f = ~ f_hypergeometric(gene_list1 = df_modules_filt %>% dplyr::filter(module == .x) %>% pull(ensembl_gene_id), 
                            gene_list2 = my_degs,
                            gene_universe = ensembl_ids
    )
)
names(l_my_degs_hypergeometric) <- sig_modules

map(l_my_degs_hypergeometric, ~ .x$overlap_genes %>% length)
map(l_my_degs_hypergeometric, ~ .x$p_value)



# Akula DEGs --------------------------------------------------------------

### HYPERGEOMETRIC (overall enrichment) ###

# Pull DEG IDs
akula_degs <- df_akula_res %>% dplyr::filter(level == "gene" & comparison == "schizo_ctrl" & pvalue < 0.05) %>% pull(id)
length(akula_degs) # n = 1373
akula_degs %in% ensembl_ids %>% sum # n = 1263 in our universe

# Run hypergeometric test for each module
l_akula_degs_hypergeometric <- map(
    .x = sig_modules,
    .f = ~ f_hypergeometric(gene_list1 = df_modules_filt %>% dplyr::filter(module == .x) %>% pull(ensembl_gene_id), 
                            gene_list2 = akula_degs,
                            gene_universe = ensembl_ids
    )
)
names(l_akula_degs_hypergeometric) <- sig_modules

map(l_akula_degs_hypergeometric, ~ .x$overlap_genes %>% length)
map(l_akula_degs_hypergeometric, ~ .x$p_value)


# Consensus DEGs ----------------------------------------------------------

### HYPERGEOMETRIC (overall enrichment) ###

# Pull consensus DEG symbols
consensus_degs <- df_consensus_degs %>% dplyr::filter(!is.na(ensembl_gene_id)) %>% pull(ensembl_gene_id)
length(consensus_degs) # n = 161 (154 with ensembl IDs)
consensus_degs %in% ensembl_ids %>% sum # n = 154 in our universe

# Run hypergeometric test for each module
l_consensus_degs_hypergeometric <- map(
    .x = sig_modules,
    .f = ~ f_hypergeometric(gene_list1 = df_modules_filt %>% dplyr::filter(module == .x) %>% pull(ensembl_gene_id), 
                            gene_list2 = consensus_degs,
                            gene_universe = ensembl_ids
    )
)
names(l_consensus_degs_hypergeometric) <- sig_modules

map(l_consensus_degs_hypergeometric, ~ .x$overlap_genes %>% length)
map(l_consensus_degs_hypergeometric, ~ .x$p_value)


# Common variants ---------------------------------------------------------

### HYPERGEOMETRIC (overall enrichment) ###

# Pull common variant IDs
common_variants <- df_scz_gwas %>% pull(gene_ensembl) %>% unique
length(common_variants) # n = 629
common_variants %in% ensembl_ids %>% sum # n = 426 in our universe

# Run hypergeometric test for each module
l_common_variants_hypergeometric <- map(
    .x = sig_modules,
    .f = ~ f_hypergeometric(gene_list1 = df_modules_filt %>% dplyr::filter(module == .x) %>% pull(ensembl_gene_id), 
                            gene_list2 = common_variants,
                            gene_universe = ensembl_ids
    )
)
names(l_common_variants_hypergeometric) <- sig_modules

map(l_common_variants_hypergeometric, ~ .x$overlap_genes %>% length)
map(l_common_variants_hypergeometric, ~ .x$p_value)


# Rare variants ---------------------------------------------------------

### HYPERGEOMETRIC (overall enrichment) ###

# Pull common variant IDs
rare_variants <- df_exome_urvs %>% pull(ensembl_gene_id) %>% unique
length(rare_variants) # n = 10
rare_variants %in% ensembl_ids %>% sum # n = 9 in our universe

# Run hypergeometric test for each module
l_rare_variants_hypergeometric <- map(
    .x = sig_modules,
    .f = ~ f_hypergeometric(gene_list1 = df_modules_filt %>% dplyr::filter(module == .x) %>% pull(ensembl_gene_id), 
                            gene_list2 = rare_variants,
                            gene_universe = ensembl_ids
    )
)
names(l_rare_variants_hypergeometric) <- sig_modules

map(l_rare_variants_hypergeometric, ~ .x$overlap_genes %>% length)
map(l_rare_variants_hypergeometric, ~ .x$p_value)



# Post-hoc: What modules are common and rare variants overrepresented in? -----------

## Curate vectors of gene names to test into a list
gene_lists <- list(
    "current_degs" = my_degs,
    "akula_degs" = akula_degs,
    "consensus_degs" = consensus_degs,
    "common_variants" = common_variants,
    "rare_variants" = rare_variants
)


## Run a hypergeometric test on each combination of benchmark list and WGCNA module
## to identify which WGCNA modules the gene lists are overrepresented in
df_gene_list_module_hypergeometric <- expand_grid(
    gene_list = names(gene_lists),
    module = names(module_colors)
) %>% 
    mutate(
        hypergeometric_res = map2(
            .x = gene_list,
            .y = module,
            .f = ~ f_module_hypergeometric(
                gene_list = gene_lists[[.x]], 
                gene_list_name = .x, 
                wgcna_module = .y, 
                gene_universe = ensembl_ids
            ) %>% 
                dplyr::select(-benchmark_gene_list, -module)
        )
    ) %>% 
    unnest(cols = c(hypergeometric_res)) %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>% 
    arrange(log10(p_adj)) %>% 
    mutate(module = factor(module, levels = names(module_colors))) %>% 
    
    # specify modules to label
    mutate(module_label = ifelse(p_adj < 0.05, paste0("n = ", overlap_n), ""))

## Plot results

l_gene_list_module_hypergeometric_plots <- map(
    .x = names(gene_lists)[3:5],
    .f = ~ f_plot_module_hypergeometric(df_gene_list_module_hypergeometric %>% dplyr::filter(gene_list == .x),
                                       y_axis = if (str_detect(.x, "degs")) {"p-adj"} else {"p-value"},
                                       plot_title = .x %>% str_replace("_", " "),
                                       module_text_size = 10,
                                       subset_x_labs = TRUE,
                                       include_module_labels = TRUE,
                                       legend_position = "none",
                                       p_value_threshold = 0.05,
                                       include_threshold_label = TRUE
    )
    
)
wrap_plots(l_gene_list_module_hypergeometric_plots, nrow = 1) &
    plot_annotation(title = "Module overrepresentation of benchmark gene lists")

## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "WGCNA_gene_list_hypergeometric", .x), width = 12, height = 4)
)

