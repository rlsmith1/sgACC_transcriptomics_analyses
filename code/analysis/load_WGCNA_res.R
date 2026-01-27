
#==============================================================================#
# Load gene-level module assignments
#==============================================================================#

## Note: This script can only be run after 02_WGCNA


base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/load_genes.R"))

## Set parameters for module set of interest
soft_power <- 3
minimum_size <- 40
tree_cut_height <- 0.98

## Load WGCNA results
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

## Filter full module assignment tibble for set of interest
df_modules_filt <- df_modules %>% 
    dplyr::filter(min_size == minimum_size & cut_height == tree_cut_height) %>% 
    tidyr::unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module) %>% 
    
    # rename modules to specify gene-level
    mutate(module = paste0("geneM", module) %>% factor(levels = paste0("geneM", levels(module)))
    )

## Define colors for plotting
module_colors <- df_modules_filt %>% dplyr::select(module, color) %>% distinct %>% deframe
