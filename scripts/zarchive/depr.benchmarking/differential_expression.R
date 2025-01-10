
########################################################################################

# Benchmark the functional enrichment of Nirmala's DE results (genes & transcripts)

########################################################################################


### SETUP ### ------

## Libraries
library(tidyverse)
library(cowplot)
library(clusterProfiler)

## Set project directory
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"

## Set figures directory
figures_dir <- paste0(base_dir, "outputs/benchmarking_res/")

## Data
load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res

# genes
prefix <- "08Mar2024_GENES_qSVAgeSexRaceGC"
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress

## Plot theme
theme_set(
    theme_cowplot() +
        theme(plot.title = element_text(size = 12),
              strip.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 8),
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 8),
              legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margin
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
        )
)


### HIGHLIGHT DEGS ### ------

df_akula_res %>% 
    filter(comparison == "schizo_ctrl" & level == "gene") %>% 
    
    filter(padj < 0.05)

### GO FUNCTIONAL ENRICHMENT OF GENES ### ------

## Set gene universe using our data
gene_universe <- df_vsd_regress %>% 
    dplyr::select(-sample) %>% 
    colnames

## Run functional enrichment
deg_go <- enrichGO(gene = df_akula_res %>% 
             filter(comparison == "schizo_ctrl" & level == "gene" & pvalue < 0.05) %>% 
             pull(id),
         OrgDb = "org.Hs.eg.db",
         universe = gene_universe,
         keyType = "ENSEMBL",
         ont = "ALL",
         pAdjustMethod = "BH",
         pvalueCutoff = 1,
         qvalueCutoff = 1,
         readable = TRUE
)

## Plot results
deg_go %>% 
    as_tibble()

### GSEA FUNCTIONAL ENRICHMENT OF GENES ### ------


### OVERLAP WITH PUBLISHED SCZ DEGS ### ------

### OVERLAP WITH PUBLISHED SCZ GWAS HITS ### ------

### OVERLAP WITH PUBLISHED CELL TYPE MARKERS ### ------


