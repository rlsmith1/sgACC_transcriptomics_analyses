
########################################################################################

# Characterize transcript-level GRCCA results with risk genes & functional enrichments

########################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/Fig5_GRCCAtranscriptsRes/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WTCNA_res.R"))


## Load GRCCA data
load(paste0(base_dir, "objects/GRCCA_transcript_results.Rdata")) # df_lvs, df_y_res, df_x_res


### Get ranked gene list ###
my_grcca_ranked <- df_x_res_transcripts %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(pearsons_r = median(pearsons_r)) %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    deframe


### Cell-type GSEA ###
df_grcca_fgsea_cell <- fgsea(pathways = cell_types, 
                             stats = my_grcca_ranked, 
                             eps = 0
)  %>% 
    as_tibble() %>% 
    dplyr::rename("cell_type" = "pathway") %>% 
    arrange(-abs(NES))
df_grcca_fgsea_cell %>% 
    mutate(sig = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "")) %>% 
    
    ggplot(aes(x = NES, y = reorder(cell_type, NES))) +
    geom_col(aes(fill = NES), color = "gray", shape = 21) +
    geom_text(aes(label = sig, x = NES), vjust = 0.75, size = 4, check_overlap = TRUE) +
    scale_fill_gradientn(colors = rev(brewer.rdbu(100)), limits = c(-3.39, 3.39)) +
    coord_cartesian(clip = "off") +
    labs(x = "NES", y = NULL,
         title = "B | Cell-type enrichment") +
    theme(legend.position = "none")


### GO GSEA ###
df_grcca_fgsea_go <- fgsea(pathways = go_pathways, 
                           stats = my_grcca_ranked, 
                           eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("go_id" = "pathway") %>% 
    arrange(-abs(NES)) %>% 
    left_join(df_go_terms) %>% # add names of GO terms
    dplyr::select(go_id, term, ontology, everything())


### Benchmark lists ###
df_grcca_fgsea_benchmark <- fgsea(pathways = benchmarking_lists, 
                                  stats = my_grcca_ranked, 
                                  eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("benchmark_list" = "pathway") %>% 
    arrange(-abs(NES))


### Export enrichment results as table ###
write_xlsx(
    list(
        "A.Cell-type" = df_grcca_fgsea_cell %>% 
            mutate(genes = map_chr(leadingEdge, ~ f_collapse_genes(.x, index = 30000))) %>% 
            dplyr::select(-leadingEdge) %>% unnest(cols = c(genes)),
        "B.GO" = df_grcca_fgsea_go %>% 
            dplyr::filter(padj < 0.05) %>%
            mutate(genes = map_chr(leadingEdge, ~ f_collapse_genes(.x, index = 30000))) %>% 
            dplyr::select(-leadingEdge) %>% unnest(cols = c(genes)),
        "C.Risk genes" = df_grcca_fgsea_benchmark %>% 
            mutate(genes = map_chr(leadingEdge, ~ f_collapse_genes(.x, index = 30000))) %>% 
            dplyr::select(-leadingEdge) %>% unnest(cols = c(genes))
    ), 
    path = paste0(tables_dir, "SXX.GRCCA_transcripts_GSEA_enrichments.xlsx")
) # exclude for now?

### Save objects for figures ###
save(df_grcca_fgsea_cell, df_grcca_fgsea_go, df_grcca_fgsea_benchmark,
     file = paste0(analysis_objects_dir, "GRCCA_transcripts_enrichment_benchmark_res.Rdata")
)

