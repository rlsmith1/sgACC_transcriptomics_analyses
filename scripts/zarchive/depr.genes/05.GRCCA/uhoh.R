

df_rcca_x_res <- read_xlsx(paste0(tables_dir, "gene_RCCA_results.xlsx"), sheet = 3)
df_grcca_x_res <- read_xlsx(paste0(tables_dir, "gene_GRCCA_results.xlsx"), sheet = 3)

df_rcca_grcca_res <- df_grcca_x_res %>% 
    dplyr::select(ensembl_gene_id, gene_symbol, z_score, pearsons_r, p_adj) %>% 
    dplyr::rename_at(3:5, ~ paste0("grcca_", .)) %>% 
    left_join(
        df_rcca_x_res %>% 
            dplyr::select(ensembl_gene_id, gene_symbol, z_score, pearsons_r, p_adj) %>% 
            dplyr::rename_at(3:5, ~ paste0("rcca_", .)),
        by = join_by(ensembl_gene_id, gene_symbol)
        
    )
df_rcca_grcca_res %>% 
    ggplot(aes(x = rcca_pearsons_r, y = grcca_pearsons_r)) +
    geom_point(aes(color = rcca_pearsons_r), size = 0.75) +
    geom_abline(lty = 2, color = "black") +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_regline_equation(label.y.npc = "top", size = 3) +
    stat_cor(label.y.npc = "top", size = 3, vjust = 3) +
    scale_color_gradientn(colors = brewer.rdbu(100) %>% rev, guide = "none") +
    labs(x = "RCCA structure correlation", y = "GRCCA structure correlation",
         title = "Relationship between RCCA and GRCCA \ngene weight structure correlations")


