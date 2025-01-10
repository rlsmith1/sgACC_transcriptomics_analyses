
######################################################################

# Gene ontology and gene set enrichment analysis of current study DEGS

######################################################################

### Set directories & load data ###

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/Fig2_DGEres/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_genes.R"))
load(paste0(base_dir, "objects/DE_results.RDS")) # df_de_res (generated in A.DGE_analysis.R)


### Get ranked gene list ###
my_degs_ranked <- df_de_res %>% 
    arrange(-stat) %>% 
    dplyr::select(ensembl_gene_id, log2fold_change) %>% 
    deframe


### Cell-type GSEA ###
df_de_fgsea_cell <- fgsea(pathways = cell_types, 
                               stats = my_degs_ranked, 
                               eps = 0
)  %>% 
    as_tibble() %>% 
    dplyr::rename("cell_type" = "pathway") %>% 
    arrange(-abs(NES))


### GO GSEA ###
df_de_fgsea_go <- fgsea(pathways = go_pathways, 
                        stats = my_degs_ranked, 
                        eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("go_id" = "pathway") %>% 
    arrange(-abs(NES)) %>% 
    left_join(df_go_terms) %>% # add names of GO terms
    dplyr::select(go_id, term, ontology, everything())


### Benchmark lists ###
df_de_fgsea_benchmark <- fgsea(pathways = benchmarking_lists, 
                               stats = my_degs_ranked, 
                               eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("benchmark_list" = "pathway") %>% 
    arrange(-abs(NES))


### Save objects for figures ###
save(df_de_fgsea_cell, df_de_fgsea_go, df_de_fgsea_benchmark,
     file = paste0(analysis_objects_dir, "DE_enrichment_benchmark_res.Rdata")
)

