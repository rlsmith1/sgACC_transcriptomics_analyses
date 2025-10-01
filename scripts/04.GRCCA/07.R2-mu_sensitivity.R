
#------------------------------------------------------------------------------#
# GENE-LEVEL: Run GRCCA across mu values to demonstrate sensitivity of the main LV
# to the group penalty (mu∈{0.001 0.01 0.5 0.9 0.999}; quasi-logarithmic scale)
#------------------------------------------------------------------------------#


## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/response_to_reviewers/round2/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
#analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WGCNA_res.R"))



# Load mu sensitivity results ---------------------------------------------

## Identify analysis directories for separate diagnosis runs
project_dir <- paste0("08Mar2024_GENES_qSVAgeSexRaceGC_sft3_minSize40_cutHeight0.98_ControlBDMDDSCZ_8MCA_regressBrainWeight_WITHGRAY/")
cca_dir <- paste0(base_dir, "RCCA_toolkit/GENES/", project_dir)

## Define variables altered for analysis
analysis_type <- "grcca"
VARx <- "0.1_1"
mus <- c(0.001, 0.01, 0.1, 0.5, 0.9, 0.999) # 0.1 is the original value


## Load GRCCA results for each
l_grcca_res <- 
    map(
        .x = mus,
        .f = ~ f_load_cca_res(
            cca_directory = cca_dir,
            level = "gene",
            analysis_type = analysis_type,
            mu = .x,
            VARx = VARx,
            include_Cmat = FALSE,
            rename_covariates = FALSE,
            all_effects = FALSE
        )
    )
names(l_grcca_res) <- paste0("mu", mu)



