
########################################################################################

# Assess functional, cell-type, and risk gene enrichments of RCCA results

########################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/supplement/")
tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WGCNA_res.R"))


## Load RCCA data
load(paste0(base_dir, "objects/RCCA_results.Rdata")) # df_lvs, df_y_res, df_x_res


### Get ranked gene list ###
my_rcca_ranked <- df_x_res_rcca %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    deframe


### Cell-type GSEA ###
df_rcca_fgsea_cell <- fgsea(pathways = cell_types, 
                          stats = my_rcca_ranked, 
                          eps = 0
)  %>% 
    as_tibble() %>% 
    dplyr::rename("cell_type" = "pathway") %>% 
    arrange(-abs(NES))


### GO GSEA ###
df_rcca_fgsea_go <- fgsea(pathways = go_pathways, 
                        stats = my_rcca_ranked, 
                        eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("go_id" = "pathway") %>% 
    arrange(-abs(NES)) %>% 
    left_join(df_go_terms) %>% # add names of GO terms
    dplyr::select(go_id, term, ontology, everything())


### Benchmark lists ###
df_rcca_fgsea_benchmark <- fgsea(pathways = benchmarking_lists, 
                               stats = my_rcca_ranked, 
                               eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("benchmark_list" = "pathway") %>% 
    arrange(-abs(NES))


### Save objects for figures ###
save(df_rcca_fgsea_cell, df_rcca_fgsea_go, df_rcca_fgsea_benchmark,
     file = paste0(analysis_objects_dir, "RCCA_enrichment_benchmark_res.Rdata")
)

