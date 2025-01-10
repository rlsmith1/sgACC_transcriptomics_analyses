
########################################################################################

# Load data and module assignments to run GRCCA & analyze outputs

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
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(ggh4x)
library(ggVennDiagram)
library(tidytext)
library(writexl)

#library(ggchicklet)


# set theme for plots -----------------------------------------------------

theme_set(theme_bw() +
              theme(plot.title = element_text(size = 11),
                    axis.title = element_text(size = 10),
                    axis.text = element_text(size = 10),
                    strip.text = element_text(size = 10),
                    legend.title = element_text(size = 8),
                    legend.text = element_text(size = 8)
              )
)

# DX COLORS
dx_colors <- c("#0072B2", "#E69F00", "#009E73", "#9966FF")
names(dx_colors) <- c("Control", "BD", "MDD", "SCZ")

# ONTOLOGY COLORS
ontology_colors <- c("BP" = "#F8766D", "CC" = "#00BA38", "MF" = "#619CFF")

# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC_CVq1"

soft_power <- 2
minimum_size <- 35
tree_cut_height <- 0.988

# LOAD OBJECTS 
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

# IDENTIFY MODULE SET OF INTEREST  
df_modules_filt <- df_modules %>% 
    filter(min_size == minimum_size & cut_height == tree_cut_height) %>% 
    tidyr::unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module) %>% 
    
    # rename modules to specify gene-level
    mutate(module = paste0("transcriptM", module) %>% factor(levels = paste0("transcriptM", levels(module)))
    )

# COVARIATES
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_clean

# DRUG MCA RESULTS
load(paste0(base_dir, "objects/drug_MCA_results.Rdata"))

# TRANSCRIPT SYMBOL TO ENSEMBL TRANSCRIPT ID MAPPING
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

# MAP TRANSCRIPTS TO GENES
df_modules_filt <- df_modules_filt %>% 
    left_join(df_gene_to_transcript) %>% 
    dplyr::select(ensembl_transcript_id, ensembl_gene_id, gene_symbol, module, color) %>% 
    distinct()

# DEFINE COLORS FOR PLOTTING
module_colors <- df_modules_filt %>% dplyr::select(module, color) %>% distinct %>% deframe

