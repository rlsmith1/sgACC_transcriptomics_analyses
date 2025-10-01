
#################################################################################################

# Load libraries, data, and functions that generalize to both gene- and transcript-level analyses

#################################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(tidymodels)
library(tidytext)

library(janitor)

library(stm)

library(patchwork)
library(ggrepel)
library(ggh4x)
library(ggpubr)
library(ggwordcloud)
library(ggvenn)
library(ggarchery)
library(ggsci)
library(ComplexUpset)
library(GGally)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(pals)

library(fgsea)
library(biomaRt)
library(DESeq2)
library(WGCNA)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(GO.db)
library(EnsDb.Hsapiens.v79) # hg38
library(ensembldb)
library(gtexr)

library(rhdf5)

library(readxl)
library(writexl)

filter <- dplyr::filter
select <- dplyr::select
count <- dplyr::count


# Load data ---------------------------------------------------------------


### CURRENT STUDY DATA ###

load(paste0(base_dir, "objects/22Jan2025_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/drug_MCA_results.Rdata")) # df_ind_loadings, df_var_loadings


### CELL TYPE DATA ###

load(paste0(base_dir, "objects/cell_type_data.RDS")) # df_cell_type

df_lake_cell_type <- df_cell_type %>% 
    # filter for Lake et al 2018 cell types
    dplyr::filter(str_detect(type, "lake")) %>% 
    mutate(type = str_remove(type, "_lake")) %>% 
    dplyr::filter(!is.na(ensembl_gene_id) & type != "synapses") %>% 
    # format cell type names for plots
    mutate(type = str_to_title(type)) %>%
    mutate(type = ifelse(type == "Opc", toupper(type), type))

cell_types <- df_lake_cell_type %>% 
    group_by(type) %>% 
    summarise(named_list = list(ensembl_gene_id)) %>% 
    deframe


### GENE ONTOLOGY ###

## List of all GO terms, each containing a vector of genes
go_pathways <- mapIds(org.Hs.eg.db, keys(org.Hs.eg.db, "GO"), "ENSEMBL", "GO", multiVals = "list")

## Full GO path annotations
df_go_terms <- AnnotationDbi::select(GO.db, names(go_pathways), c("TERM", "ONTOLOGY"), "GOID", multiVals = "first") %>% 
    as_tibble() %>% 
    clean_names() %>% 
    dplyr::rename("go_id" = "goid")



# Load benchmarking risk gene lists ---------------------------------------


## SCZ GWAS DATA (broad set)
df_scz_broad_variants <- readxl::read_excel(paste0(base_dir, "data/trubetskoy_2023_scz_gwas.xlsx"),
                                            sheet = 3
)
scz_broad_variants <- df_scz_broad_variants %>% dplyr::filter(str_detect(gene_ensembl, "^ENSG0")) %>% pull(gene_ensembl) %>% unique

# sheet 2 = 95% credible set
# sheet 3 = 95% credible set k <= 3.5 (broad fine-mapped set)
# sheet 4 = 95% credible set <= 5 SNPs
# sheet 5 = 95% credible set <= 5 PIP

## Prioritized SCZ GWAS variants
df_scz_prioritized_variants <- readxl::read_excel(paste0(base_dir, "data/Trubetskoy_prioritized_genes.xlsx"))
scz_prioritized_variants <- df_scz_prioritized_variants %>% dplyr::filter(!is.na(ensembl_gene_id)) %>% pull(ensembl_gene_id) %>% unique

## SCZ rare variants
df_scz_rare_variants <- readxl::read_excel(paste0(base_dir, "data/singh_2022_exome_URVs.xlsx"))
scz_rare_variants <- df_scz_rare_variants %>% pull(ensembl_gene_id) %>% unique

## Bipolar significant MAGMA genes
df_bd_gwas <- read_xlsx(paste0(base_dir, "data/mullins_2021_BD_GWAS.xlsx"), sheet = "Table S4") %>% 
    row_to_names(1) %>% 
    clean_names %>% 
    filter(!is.na(ensembl_id))
bd_common_variants <- df_bd_gwas %>% pull(ensembl_id) %>% unique

## MDD MAGMA genes
df_mdd_gwas <- read_xlsx(paste0(base_dir, "data/howard_2019_MDD_GWAS.xlsx")) %>% 
    row_to_names(1) %>% 
    clean_names %>% 
    mutate(ensembl_gene_id = mapIds(org.Hs.eg.db, keys = gene_name, "ENSEMBL", "SYMBOL", multiVals = "list"))
mdd_common_variants <- df_mdd_gwas %>% pull(ensembl_gene_id) %>% unlist %>% unique
mdd_common_variants <- mdd_common_variants[!is.na(mdd_common_variants)]

## ASD MAGMA genes
df_asd_gwas <- read_xlsx(paste0(base_dir, "data/matoba_2020_ASD_GWAS.xlsx")) %>% 
    clean_names
asd_common_variants <- df_asd_gwas %>% pull(gene) %>% unique


## Combine benchmarking lists into list
benchmarking_lists <- list(
    "SCZ (common, broad)" = scz_broad_variants,
    "SCZ (common, prioritized)" = scz_prioritized_variants,
    "SCZ (rare)" = scz_rare_variants,
    "BD" = bd_common_variants,
    "MDD" = mdd_common_variants,
    "ASD" = asd_common_variants
)



# Group covariates by type ------------------------------------------------


technical_covariates <- c(
    "mapped_percent", "gc_percent", "five_prime_three_prime_bias", "rin_acsg", "rna_extraction_batch",
    "library_batch", "pmi", "pmi_confidence", "source", "ph", "max_rin", "max_rine"
)
biological_covariates <- c(
    "age_death", "sex_at_birth", "race", "bmi", "height", "weight",
    "marital_status", "manner_death", "education", "suicide", "brain_weight"
)
drug_covariates <- c(
    "smoker", "nicotine_cotinine", "alcohol", "opioids", "cannabinoids", "beta_blockers",
    "major_stimulants_cocaine_included", "minor_stimulants", "anticholinergics", "antidepressants",
    "anti_epileptics", "anti_histamines", "antipsychotics", "mood_stabilizers", "sedative_hypnotic_anxiolitics", "benzos",
    "other_psychotropic_drug", "non_psychiatric"
)


# Source functions --------------------------------------------------------

## for plots
source(paste0(base_dir, "functions/plot_functions.R"))

## genes hypergeometric
source(paste0(base_dir, "functions/hypergeometric_functions.R"))

## load CCA res
source(paste0(base_dir, "functions/f_load_CCA_res.R"))

## Helper functions

# Collapse string of genes from GSEA results
f_collapse_genes <- function(x, index = 10000) {
    x <- paste(x, collapse = "|")
    x <- paste0(substr(x, 1, index), "...")
}



# Set plot theme ----------------------------------------------------------

## Plot aesthetics
theme_set(theme_cowplot() +
              theme(plot.title = element_text(size = 11),
                    axis.title = element_text(size = 10),
                    axis.text = element_text(size = 9),
                    strip.text = element_text(size = 10),
                    legend.title = element_text(size = 9),
                    legend.text = element_text(size = 8),
                    plot.margin = margin(t = 0, r = 5, l = 0, b = 0)
              )
)

## Figure dimensions
fig_width <- 2.5
fig_height <- 2

## Dx colors
dx_colors <- c(
    "Control" = "#0072B2", 
    "BD" = "#E69F00", 
    "MDD" = "#009E73", 
    "SCZ" = "#9966FF"
)

