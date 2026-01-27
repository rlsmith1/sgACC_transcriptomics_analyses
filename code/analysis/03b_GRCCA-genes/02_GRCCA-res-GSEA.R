
#==============================================================================#
# Characterize the biological enrichments of GRCCA data
#==============================================================================#

## Run GSEA on GRCCA results
## Note that the input data (cell types, GO, and risk gene benchmark lists)
## are loaded & collated in ~/setup.R


#----Set directories & load data----#

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(scripts_dir, "load_WGCNA_res.R"))
load(paste0(objects_dir, "GRCCA_results.Rdata")) # df_lvs, df_y_res, df_x_res


#----Get ranked gene list for GSEA----#

my_grcca_ranked <- df_x_res %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    deframe


#----Cell type----#

df_grcca_fgsea_cell <- fgsea(pathways = cell_types, 
                          stats = my_grcca_ranked, 
                          eps = 0
)  %>% 
    as_tibble() %>% 
    dplyr::rename("cell_type" = "pathway") %>% 
    arrange(-abs(NES))


#----GO----#

df_grcca_fgsea_go <- fgsea(pathways = go_pathways, 
                        stats = my_grcca_ranked, 
                        eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("go_id" = "pathway") %>% 
    arrange(-abs(NES)) %>% 
    left_join(df_go_terms) %>% # add names of GO terms
    dplyr::select(go_id, term, ontology, everything())


#----Risk genes----#

df_grcca_fgsea_benchmark <- fgsea(pathways = benchmarking_lists, 
                               stats = my_grcca_ranked, 
                               eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("benchmark_list" = "pathway") %>% 
    arrange(-abs(NES))


#----Save----#

save(df_grcca_fgsea_cell, df_grcca_fgsea_go, df_grcca_fgsea_benchmark,
     file = paste0(objects_dir, "GRCCA_enrichment_benchmark_res.Rdata")
)

