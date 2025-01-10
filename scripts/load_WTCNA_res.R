
### NOTE: This script can only be run after scripts/06.GRCCA_transcripts/00.case_control_WTCNA.R is run! ###

################################################################################

# Load transcript-level module assignments

################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/load_transcripts.R"))

## Define parameters for module set of interest
soft_power <- 2
minimum_size <- 35
tree_cut_height <- 0.988

## Load module data
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

## Filter module data for correct module set
df_modules_filt <- df_modules %>% 
    dplyr::filter(min_size == minimum_size & cut_height == tree_cut_height) %>% 
    tidyr::unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module) %>% 
    
    # rename modules to specify gene-level
    mutate(module = paste0("transcriptM", module) %>% factor(levels = paste0("transcriptM", levels(module)))
    )

## Define module colors for plot
module_colors <- df_modules_filt %>% dplyr::select(module, color) %>% distinct %>% deframe
