#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load transcript-level GRCCA results (with & without Cmat, or Dim2 regression)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/Fig5_GRCCAtranscriptsRes/")
tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/Fig5_GRCCAtranscriptsRes/")

source(paste0(base_dir, "scripts/full_analysis_scripts/scripts/load_WTCNA_res.R"))


# Load data ---------------------------------------------------------------

# 08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC_CVq1_sft2_minSize35_cutHeight0.988_8MCA_ControlBDMDDSCZ_regressDim2_WITHGray

## Identify CCA directory
dx_groups_to_analyze <- c("Control", "BD", "MDD", "SCZ") # select diagnostic groups of interest
n_mca_dim <- 8
type <- "OVERALLsubs/"

## Read in GRCCA analysis results
analysis_type <- "grcca"
VARx <- "0.1_1"

# with Cmat --> DOES NOT include Dim 2 drugs in Y results! (regress!)
project_dir <- paste0(
    prefix,
    "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "_",
    paste0(dx_groups_to_analyze, collapse = ""), "_",
    n_mca_dim, "MCA_regressDim2_WITHGRAY/"
)
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)
l_grcca_res_withCmat <- f_load_cca_res(
    cca_directory = cca_dir,
    level = "transcript",
    analysis_type = analysis_type,
    mu = 0.1,
    VARx = VARx,
    include_Cmat = TRUE
)

# no Cmat --> includes Dim 2 drugs in Y results!
project_dir <- paste0(
    prefix,
    "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "_",
    paste0(dx_groups_to_analyze, collapse = ""), "_",
    n_mca_dim, "MCA_regressBrainWeight_WITHGRAY/"
)
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)
l_grcca_res_noCmat <- f_load_cca_res(
    cca_directory = cca_dir,
    level = "transcript",
    analysis_type = analysis_type,
    mu = 0.1,
    VARx = VARx,
    include_Cmat = FALSE
)



# Assign results to tibbles -----------------------------------------------


## Assign results to tibbles; make sure SCZ is positive for interpretability
df_results <- bind_rows(
    l_grcca_res_withCmat$model_results %>% mutate(model = "regressDim2", .before = 1),
    l_grcca_res_noCmat$model_results %>% mutate(model = "includeDim2", .before = 1)
)
df_lvs <- bind_rows(
    l_grcca_res_withCmat$latent_variables %>% mutate(model = "regressDim2", .before = 1),
    l_grcca_res_noCmat$latent_variables %>% mutate(model = "includeDim2", .before = 1)
)
df_y_res <-  bind_rows(
    l_grcca_res_withCmat$y_res %>% mutate(model = "regressDim2", .before = 1),
    l_grcca_res_noCmat$y_res %>% mutate(model = "includeDim2", .before = 1) 
)
df_x_res <- bind_rows(
    l_grcca_res_withCmat$x_res %>% mutate(model = "regressDim2", .before = 1),
    l_grcca_res_noCmat$x_res %>% mutate(model = "includeDim2", .before = 1) 
) %>% 
    
    left_join(df_transcript_to_gene, by = join_by(ensembl_transcript_id)) %>% 
    left_join(df_modules_filt %>% dplyr::select(-c(mod_set, color)), by = join_by(ensembl_transcript_id)) %>% 
    dplyr::select(model, ensembl_transcript_id, transcript_symbol, module, everything()) %>% 
    mutate(transcript_symbol = ifelse(is.na(transcript_symbol), ensembl_transcript_id, transcript_symbol)) %>% 
    distinct()



# Run benchmarking tests on both sets of results --------------------------

### Get ranked gene lists ###
regressDim2_grcca_ranked <- df_x_res %>% 
    filter(model == "regressDim2") %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(pearsons_r = median(pearsons_r)) %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    deframe
includeDim2_grcca_ranked <- df_x_res %>% 
    filter(model == "includeDim2") %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(pearsons_r = median(pearsons_r)) %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    deframe
ranked_grcca_lists <- list("regressDim2" = regressDim2_grcca_ranked, "includeDim2" = includeDim2_grcca_ranked)


### Cell-type GSEA ###
df_grcca_fgsea_cell <- map_dfr(
    .x = c("regressDim2", "includeDim2"),
    .f = ~ fgsea(pathways = cell_types, 
                 stats = ranked_grcca_lists[[.x]], 
                 eps = 0
    )  %>% 
        as_tibble() %>% 
        dplyr::rename("cell_type" = "pathway") %>% 
        arrange(-abs(NES)) %>% 
        mutate(model = .x, .before = 1)
)

### GO GSEA ###
df_grcca_fgsea_go <- map_dfr(
    .x = c("regressDim2", "includeDim2"),
    .f = ~ fgsea(pathways = go_pathways, 
                 stats = ranked_grcca_lists[[.x]], 
                 eps = 0
    )  %>% 
        as_tibble() %>% 
        dplyr::rename("go_id" = "pathway") %>% 
        arrange(-abs(NES)) %>% 
        left_join(df_go_terms) %>% # add names of GO terms
        dplyr::select(go_id, term, ontology, everything()) %>% 
        mutate(model = .x, .before = 1)
)

### Benchmark lists ###
df_grcca_fgsea_benchmark <- map_dfr(
    .x = c("regressDim2", "includeDim2"),
    .f = ~ fgsea(pathways = benchmarking_lists, 
                 stats = ranked_grcca_lists[[.x]], 
                 eps = 0
    )  %>% 
        as_tibble() %>% 
        dplyr::rename("benchmark_list" = "pathway") %>% 
        arrange(-abs(NES)) %>% 
        mutate(model = .x, .before = 1)
)


# SAVE for figures --------------------------------------------------------


save(
    df_results, df_lvs, df_y_res, df_x_res,
    file = paste0(analysis_objects_dir, "transcript_GRCCA_res.Rdata")
)

save(
    df_grcca_fgsea_cell, df_grcca_fgsea_go, df_grcca_fgsea_benchmark,
    file = paste0(analysis_objects_dir, "GRCCA_transcripts_enrichment_benchmark_res.Rdata")
)

