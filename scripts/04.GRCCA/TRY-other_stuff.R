
################################################################################

# GENE-LEVEL: Other things to enhance GRCCA results...

################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/updated_figures/figures/Fig4_characterizeGRCCA/")
tables_dir <- paste0(base_dir, "outputs/updated_figures/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/updated_figures/objects/Fig4_characterizeGRCCA/")

source(paste0(base_dir, "scripts/full_analysis_scripts/scripts/load_WGCNA_res.R"))

## Load GRCCA data
load(paste0(base_dir, "objects/GRCCA_results.Rdata")) # df_lvs, df_y_res, df_x_res



# gene-eQTL mapping -------------------------------------------------------


df_eqtls <- read.table("~/Downloads/DER-08a_hg38_eQTL.significant.txt") %>% 
    row_to_names(1) %>%
    as_tibble()


df_eqtls_grcca <- df_eqtls %>% 
    mutate(ensembl_gene_id = str_remove(gene_id, "\\..*")) %>% 
    dplyr::select(ensembl_gene_id, SNP_id, nominal_pval, FDR, regression_slope, top_SNP) %>% 
    left_join(df_x_res) %>% 
    dplyr::filter(!is.na(pearsons_r))

## Summarise genes by taking median regression slope across eQTLs
df_eqtls_grcca_avg <- df_eqtls_grcca %>%
    group_by(ensembl_gene_id) %>% 
    mutate(regression_slope = as.numeric(regression_slope)) %>% 
    summarise(
        regression_slope = median(regression_slope),
        pearsons_r = median(pearsons_r)
    )

## Summarise genes by taking only the top SNP
df_eqtls_grcca_top <- df_eqtls_grcca %>%
    dplyr::filter(top_SNP == 1) %>% 
    mutate(regression_slope = as.numeric(regression_slope))

## Plot
df_eqtls_grcca_top %>% 
    #mutate(eqtl_decile = ntile(-regression_slope, 10) %>% factor) %>% 
    
    #dplyr::filter(regression_slope > 0) %>% 
    ggplot(aes(x = regression_slope, y = pearsons_r)) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_cor() #+
    facet_wrap(vars(eqtl_decile), scales = "free")



# How does residual of plot 3E relate to risk gene-ness -------------------
    
   
## Combine DE & GRCCA 
df_de_grcca_res <- df_x_res %>% 
        left_join(df_de_res)    
   

# Find equation for line of best fit
    lm_grcca_de <- lm(scale(pearsons_r) ~ scale(stat), data = df_de_grcca_res)
    
    summary(lm_grcca_de)
    
    
# Function to calculate distance from point to line
    f_dist_point_to_line <- function(m, b, x0, y0) {
        
        # Formula to calculate the distance from a point to a line y = mx + b
        distance <- (m * x0 - y0 + b) / sqrt(m^2 + 1)
        
        return(distance)
    }
    

## Calculate distance for each point to line of best fit    
    df_de_grcca_res_dist <- df_de_grcca_res %>% 
        mutate(distance = f_dist_point_to_line(
            m = coefficients(lm_grcca_de)[2],
            b = coefficients(lm_grcca_de)[1],
            x0 = scale(stat)[,1],
            y0 = scale(pearsons_r)[,1]
        )
        )
    
    df_de_grcca_res_dist %>% 
        arrange(-distance)
    
## Calculate distance for each point to identity line    
    df_de_grcca_res_dist <- df_de_grcca_res %>% 
        mutate(distance = f_dist_point_to_line(
            m = 1,
            b = 0,
            x0 = scale(stat)[,1],
            y0 = scale(pearsons_r)[,1]
        )
        )
    
    df_de_grcca_res_dist %>% 
        arrange(-distance)
    df_de_grcca_res_dist %>% 
        arrange(distance)
    
## Plot     
    df_de_grcca_res_dist %>% 
        mutate(label = ifelse(abs(distance) > 2, gene_symbol, "")) %>% 
        
        ggplot(aes(x = scale(stat), y = scale(pearsons_r))) +
        geom_hline(yintercept = 0, color = "gray") +
        geom_vline(xintercept = 0, color = "gray") +
        geom_point(aes(color = pearsons_r), size = 0.5) +
        geom_abline(lty = 2) +
        geom_smooth(method = "lm", color = "black", linewidth = 0.5) +
        geom_text_repel(aes(label = label), min.segment.length = 0) +
        stat_cor(aes(label = after_stat(r.label)), label.sep = "\n", size = 3,
                 label.y.npc = "top", vjust = 0) +
        scale_color_gradientn(colors = gene_weight_color_scale, limits = c(-0.57, 0.57))
    
## GSEA enrichments for risk genes?
    gene_distance_ranked <- df_de_grcca_res_dist %>% 
        arrange(-distance) %>% 
        dplyr::select(ensembl_gene_id, distance) %>%
        deframe
    
   fgsea(pathways = benchmarking_lists, 
                                   stats = gene_distance_ranked, 
                                   eps = 0
    ) %>% 
        as_tibble() %>% 
        dplyr::rename("benchmark_list" = "pathway") %>% 
        arrange(-abs(NES))
    

# Correlate enrichments ---------------------------------------------------

## Load enrichment data objects
   directories <- c(
       paste0(base_dir, "outputs/updated_figures/objects/Fig2_DGEres/"),
       paste0(base_dir, "outputs/updated_figures/objects/Fig4_characterizeGRCCA/")
   )
   
for (dir in directories) {
    objects <- list.files(dir)
    objects <- objects[!str_detect(objects, "null|archive")]
    for (obj in objects){
        print(paste0("loading ", obj, "...."))
        load(paste0(dir, obj))
    }
}
   
   
# Scatter plot comparison
   
   ## Cell type
   df_de_fgsea_cell %>%
       mutate(analysis = "DGE", .before = 1) %>% 
       bind_rows(
           df_grcca_fgsea_cell %>%
               mutate(analysis = "GRCCA", .before = 1)
       ) %>% 
       
       pivot_wider(id_cols = cell_type, names_from = analysis, values_from = NES) %>% 
       ggplot(aes(x = DGE, y = GRCCA)) +
       geom_point() +
       geom_text_repel(aes(label = cell_type)) +
       geom_hline(yintercept = 0, color = "gray") +
       geom_vline(xintercept = 0, color = "gray") +
       geom_smooth(method = "lm", se = FALSE, color = "black")
   
   ## Risk gene
   df_de_fgsea_benchmark %>%
       mutate(analysis = "DGE", .before = 1) %>% 
       bind_rows(
           df_grcca_fgsea_benchmark %>%
               mutate(analysis = "GRCCA", .before = 1)
       ) %>% 
       
       pivot_wider(id_cols = benchmark_list, names_from = analysis, values_from = NES) %>% 
       ggplot(aes(x = DGE, y = GRCCA)) +
       geom_point() +
       geom_text_repel(aes(label = benchmark_list)) +
       # geom_hline(yintercept = 0, color = "gray") +
       # geom_vline(xintercept = 0, color = "gray") +
       geom_smooth(method = "lm", se = FALSE, color = "black")
   
   ## GO term
   df_de_fgsea_go %>%
       mutate(analysis = "DGE", .before = 1) %>% 
       bind_rows(
           df_grcca_fgsea_go %>%
               mutate(analysis = "GRCCA", .before = 1) %>% 
               dplyr::filter(pval < 0.05)
       ) %>%
       pivot_wider(id_cols = c(term, ontology), names_from = analysis, values_from = NES) %>% 
       #mutate(label = ifelse(GRCCA < -2.75, term, "")) %>% 
       
       ggplot(aes(x = DGE, y = GRCCA)) +
       geom_point(aes(color = ontology)) +
       stat_cor() +
       geom_text_repel(aes(label = term), max.overlaps = 20) +
       geom_hline(yintercept = 0, color = "gray") +
       geom_vline(xintercept = 0, color = "gray") +
       geom_abline(lty = 2, color = "darkgray") +
       geom_smooth(method = "lm", se = FALSE, color = "black")
   
   ## Opposite sign?
   df_de_fgsea_go %>%
       mutate(analysis = "DGE", .before = 1) %>% 
       bind_rows(
           df_grcca_fgsea_go %>%
               mutate(analysis = "GRCCA", .before = 1)
       ) %>% 
       pivot_wider(id_cols = c(term, ontology), names_from = analysis, values_from = NES) %>% 
       mutate(dist = DGE - GRCCA) %>% 
       arrange(-dist) %>% 
       dplyr::filter(DGE > 0 & GRCCA < 0)
   
   df_de_fgsea_go %>%
       mutate(analysis = "DGE", .before = 1) %>% 
       bind_rows(
           df_grcca_fgsea_go %>%
               mutate(analysis = "GRCCA", .before = 1)
       ) %>% 
       pivot_wider(id_cols = c(term, ontology), names_from = analysis, values_from = NES) %>% 
       mutate(dist = DGE - GRCCA) %>% 
       arrange(dist) %>% 
       dplyr::filter(DGE < 0 & GRCCA > 0)
   
   
   
