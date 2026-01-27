
#==============================================================================#
# Gene ontology and gene set enrichment analysis of current study DEGS
#==============================================================================#

## Run GSEA on limma DGE results
## Note that the input data (cell types, GO, and risk gene benchmark lists)
## are loaded & collated in ~/setup.R


#----Set directories & load data----#

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(scripts_dir, "load_genes.R"))
load(paste0(base_dir, "objects/DE_results.Rdata")) # df_limma_res (generated in 00_DGE-analysis-limma.R)


#----Get ranked gene list for GSEA----#

my_degs_ranked <- df_limma_res %>% 
    arrange(-log_fc) %>%
    dplyr::select(ensembl_gene_id, log_fc) %>%
    # arrange(-t) %>%
    # dplyr::select(ensembl_gene_id, t) %>% 
    deframe


#----Cell type----#

df_de_fgsea_cell <- fgsea(pathways = cell_types, 
                               stats = my_degs_ranked, 
                               eps = 0
)  %>% 
    as_tibble() %>% 
    dplyr::rename("cell_type" = "pathway") %>% 
    arrange(-abs(NES))


#----GO----#

df_de_fgsea_go <- fgsea(pathways = go_pathways, 
                        stats = my_degs_ranked, 
                        eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("go_id" = "pathway") %>% 
    arrange(-abs(NES)) %>% 
    left_join(df_go_terms) %>% # add names of GO terms
    dplyr::select(go_id, term, ontology, everything())


#----Risk genes----#

df_de_fgsea_benchmark <- fgsea(pathways = benchmarking_lists, 
                               stats = my_degs_ranked, 
                               eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("benchmark_list" = "pathway") %>% 
    arrange(-abs(NES))


#----Save----#

save(df_de_fgsea_cell, df_de_fgsea_go, df_de_fgsea_benchmark,
     file = paste0(analysis_objects_dir, "DE_enrichment_benchmark_res.Rdata")
)

