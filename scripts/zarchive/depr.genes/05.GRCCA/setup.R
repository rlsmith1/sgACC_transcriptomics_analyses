
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
library(ggh4x)
library(ggpubr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(readxl)
library(writexl)
library(cowplot)


# set theme for plots -----------------------------------------------------

theme_set(theme_cowplot() +
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
prefix <- "08Mar2024_GENES_qSVAgeSexRaceGC"

soft_power <- 3
minimum_size <- 40
tree_cut_height <- 0.98

# LOAD OBJECTS 
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

# IDENTIFY MODULE SET OF INTEREST  
df_modules_filt <- df_modules %>% 
    filter(min_size == minimum_size & cut_height == tree_cut_height) %>% 
    tidyr::unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module) %>% 
    
    # rename modules to specify gene-level
    mutate(module = paste0("geneM", module) %>% factor(levels = paste0("geneM", levels(module)))
    )

module_colors <- df_modules_filt %>% dplyr::select(module, color) %>% distinct %>% deframe

