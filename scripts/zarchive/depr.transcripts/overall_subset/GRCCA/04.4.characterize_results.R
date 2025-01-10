
########################################################################################

# Characterize transcript-level GRCCA results

########################################################################################


# X-Y LV correlation ------------------------------------------------------


cor.test(df_lvs$lvx, df_lvs$lvy, method = "spearman")

df_lvs %>% 
  ggplot(aes(x = lvx, y = lvy)) +
  geom_point(aes(color = dx), alpha = 0.7, size = 3) +
  geom_smooth(method = lm, se = FALSE, color = "black") +
  stat_cor() +
  scale_color_manual(values = dx_colors) +
  guides(color = guide_legend(title = "diagnosis")) +
  labs(x = "X latent variable (x matrix • x weights)",
       y = "Y latent variable (y matrix • y weights)",
       title = "Transcript x-y latent variable correlation") +
  theme(legend.position = c(0.8, 0.15))




# Fig 5A | y weights ---------------------------------------------------------------

## WEIGHTS
df_y_res %>% 
  mutate(significant = ifelse(p_adj < 0.05 & abs(z_score) >= 2, "*", "")) %>% 

  ggplot(aes(x = reorder(covariate, weight), y = weight)) +
  geom_col(aes(fill = covariate, alpha = abs(weight)), color = "black") +
  geom_errorbar(aes(ymin = weight - sd, ymax = weight + sd), width = 0.2) +
  geom_text(aes(label = significant, y = 0), size = 10, vjust = 1.2) +
  geom_hline(aes(yintercept = 0), color = "black") +
  scale_fill_manual(values = c(dx_colors, rep("gray", 6))) +
  scale_alpha_continuous(range = c(0.1, 1)) +
  scale_x_discrete(labels = function(x) {str_replace(x, "_", " ")}) +
  ylim(c(-0.15, 0.15)) +
  labs(x = "", y = "weight", 
       title = "Covariate weights",
       caption = "Error bar shows +/- 1 sd") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 0.9))

# save to project dir
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/5A.yres_weights", .x),
                width = 4, height = 3)
)



# Fig 5B | y structure correlations -------------------------------------------


# Y STRUCTURE CORRELATIONS  
df_y_res %>% 
    mutate(significant = ifelse(p_adj < 0.05 & abs(z_score) >= 2, "*", "")) %>% 
    
    ggplot(aes(x = pearsons_r, y = reorder(covariate, pearsons_r))) +
    geom_point(aes(size = -log10(p_adj), fill = pearsons_r), shape = 21) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_text(aes(label = significant), size = 6, vjust = 0.75, hjust = -1.25, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
    scale_size_continuous(range = c(2, 6)) +
    guides(fill = guide_colorbar(title = "Pearson's r"),
           size = guide_legend(title = "-log10(p-adj)")) +
    #scale_color_manual(values = c("*" = "black", "-" = "transparent")) +
    xlim(c(-1.0, 1.0)) +
    labs(y = "", x = "Pearson's r", 
         title = "Covariate structure correlations") +
    theme(legend.position = c(0.15, 0.60), 
          legend.box = "vertical",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.30, "cm")
    )

# save to project dir
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/5B.yres_struc_cor", .x),
                width = 4, height = 4)
)



# Fig 5C | Significant X res ----------------------------------------------


### VECTOR OF ALL SIGNIFICANT GENES

df_x_res %>% 
    filter(p_adj < 0.05 & abs(z_score) >= 2) %>% # n = 257
    
    ggplot(aes(x = pearsons_r, y = reorder(ensembl_transcript_id, pearsons_r))) +
    geom_col(aes(fill = pearsons_r)) +
    geom_vline(xintercept = 0, lty = 2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
    guides(fill = guide_colorbar(title = "Pearson's r")) +
    #scale_color_manual(values = c("*" = "black", "-" = "transparent")) +
    xlim(c(-1, 1)) +
    labs(y = "Isoform", x = "Pearson's r", 
         title = "Significant isoform structure correlations",
         #title = bquote(bold("D")~" | Significant gene structure correlations"),
         caption = "n = 257"
    ) +
    theme_classic() +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          legend.position = "none"
    )

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/5C.xres_struc_cor_significant", .x), # or grcca/withGray
                  width = 4, height = 6.5)
)


### STRUCTURE CORRELATIONS OF TOP N GENES

df_x_res %>% 
    filter(p_adj < 0.05 & abs(z_score) >= 2) %>%
    arrange(-abs(pearsons_r)) %>% 
    head(25) %>% 
    mutate(direction = ifelse(pearsons_r < 0, "2neg", "1pos")) %>% 
    
    ggplot(aes(x = pearsons_r, y = reorder(transcript_symbol, pearsons_r))) +
    geom_point(aes(size = -log10(p_adj), fill = pearsons_r), shape = 21) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_wrap(vars(direction), ncol = 1, scales = "free_y") +
    force_panelsizes(rows = c(6, 19)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         limits = c(-1, 1), guide = "none") +
    scale_size_continuous(range = c(2, 6)) +
    guides(size = guide_legend(title = "-log10(p-adj)")) +
    xlim(c(-1, 1)) +
    labs(y = "", x = "Pearson's r", 
         title ="Top 25 by abs(Pearson's r)"
    ) +
    theme(strip.text = element_blank(), 
          panel.spacing = unit(3, "lines"),
          legend.position = c(0.75, 0.20), 
          legend.box = "horizontal",
          legend.key.size = unit(0.35, "cm")
    )

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/5C.xres_struc_cor_top25", .x), # or grcca/withGray
                  width = 4, height = 6.5)
)



# Fig 5D | Module hypergeometric overrepresentation ----------------------------------------------------------------


## HYPERGEOMETRIC TEST FOR SIGNIFICANT OVERLAP
transcript_universe <- df_x_res %>% 
    pull(ensembl_transcript_id)
grcca_transcripts <- df_x_res %>%
    filter(significant == 1) %>%
    pull(ensembl_transcript_id) %>% 
    unique

df_module_overlap <- tibble()
for (mod in unique(df_x_res$module))  {
    
    mod_transcripts <- df_x_res %>% 
        filter(module == mod) %>% 
        pull(ensembl_transcript_id)
    
    mod_grcca_intersect <- intersect(grcca_transcripts, mod_transcripts)
    
    p <- 1 - phyper(q = length(mod_grcca_intersect), 
                    m = length(grcca_transcripts), 
                    n = length(transcript_universe) - length(grcca_transcripts), 
                    k = length(mod_transcripts)
    )
    
    df_tmp <- tibble(module = mod,
                     overlap = length(mod_grcca_intersect),
                     mod_size = length(mod_transcripts),
                     p_value = p)
    
    df_module_overlap <- df_module_overlap %>% bind_rows(df_tmp)
    
}
df_module_overlap <- df_module_overlap %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>% 
    left_join(df_modules_filt %>% dplyr::select(module, color) %>% distinct()) %>% 
    mutate(module = factor(module, levels = unique(df_modules_filt$module))) %>% 
    mutate(p_value = ifelse(overlap <= 1, 1, p_value),
           p_adj = ifelse(overlap <= 1, 1, p_adj)
    ) %>% 
    arrange(p_value)

# PLOT
df_module_overlap %>%
    ggplot(aes(x = module, y = -log10(p_adj))) +
    geom_col(width = 0.05, fill = "black") +
    geom_point(aes(fill = module, size = overlap/mod_size), shape = 21) +
    geom_hline(aes(yintercept = -log10(0.01)), linewidth = 0.25) +
    annotate(geom = "text", label = "FDR < 0.01", x = 23, y = -log10(0.01) + 0.2) +
    scale_fill_manual(values = module_colors) +
    labs(x = "", y = "-log10(hypergeometric FDR)",
         title = "Module overrepresentation in GRCCA significant transcript list") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 0.9, size = 7)
    )

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/5D.module_hypergeometric_res", .x), # or grcca/withGray
                  width = 5, height = 4)
)



# Fig 5E | Directions of signficant genes in top modules ------------------------------------------------------------



## SHOW TRANSCRIPT DIRECTION IN EACH SIGNIFICANT MODULE
sig_modules <- df_module_overlap %>% filter(p_adj < 0.01) %>% pull(module) %>% sort()

df_x_res %>% 
    filter(module %in% sig_modules & ensembl_transcript_id %in% grcca_transcripts) %>%
    mutate(weight = ifelse(significant == 1, weight, 0)) %>% 
    arrange(module) %>% 
    mutate(ensembl_transcript_id = factor(ensembl_transcript_id, levels = .$ensembl_transcript_id)) %>% 
    
    ggplot(aes(x = ensembl_transcript_id, y = pearsons_r)) +
    geom_col(aes(fill = module, alpha = abs(pearsons_r))) +
    facet_wrap(vars(module), nrow = 1, scales = "free_x") +
    force_panelsizes(cols = c(18, 11, 59, 3)) +
    scale_fill_manual(values = module_colors) +
    scale_alpha_continuous(range = c(0.3, 1)) +
    labs(title = "Significant isoforms in overrepresented modules", 
         x = "Isoform", y = "Pearson's r") +
    guides(fill = guide_legend(nrow = 2)) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_blank(),
          legend.position = "none")

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/5E.xres_sig_mods_gene_directions", .x), # or grcca/withGray
                  width = 8, height = 4)
)



# Fig 5F | Significant module GO enrichments ------------------------------


# Load GO results
load( paste0(base_dir, "objects/", prefix,
             "_SIGNED_SFT", soft_power, 
             "_MIN_SIZE", minimum_size, 
             "_CUT_HEIGHT", tree_cut_height, "_GO_RES.RDS")
)

## PLOT
df_mods_go %>% 
    clean_names %>% 
    filter(!is.na(description) & module %in% sig_modules) %>% 
    
    # take top 5 paths by p-value per module to plot
    dplyr::group_by(module) %>% 
    arrange(module, -pvalue) %>% 
    mutate(row = row_number(),
           gene_ratio = parse(text = gene_ratio) %>% eval
    ) %>% 
    top_n(n = 5, wt = row) %>% 
    #slice_max(order_by = -weight_fisher, n = 5) %>% 
    mutate(`survives FDR` = ifelse(p_adjust < 0.05, "yes", "no")) %>% 
    mutate(description = str_wrap(description, width = 30)) %>%
    
    # plot
    ggplot(aes(x = -log10(pvalue), y = reorder_within(description, within = module, by = -pvalue))) +
    #y = reorder_within(str_wrap(description, width = 35), -pvalue, module))) +
    geom_point(aes(size = gene_ratio, fill = ontology), #color = `survives FDR`),
               color = "black", shape = 21, stroke = 1) +
    geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
    scale_size_continuous(range = c(2, 5)) +
    #scale_color_manual(values = c("yes" = "black", "no" = "transparent")) +
    scale_y_reordered() +
    facet_wrap(vars(module), scales = "free_y", ncol = 4) +
    labs(y = "", x = "log10(p-value)",
         title = "Overrepresented modules enrichments") +
    guides(size = guide_legend(title = "Gene ratio")) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title.position = "top"
    )


# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/5F.xres_sig_mods_GO", .x), # or grcca/withGray
                  width = 11, height = 3)
)


# Fig 5G | Gene set enrichment analyses ---------------------------------------------------


# GENERATE GENE VECTOR TO RUN GSEA ON (MAP TRANSCRIPT --> GENE)
ranked_gene_list <- df_x_res %>%
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(pearsons_r = median(pearsons_r)) %>% 
    arrange(-pearsons_r) %>%
    deframe 
length(ranked_gene_list) # n = 18032


# GENE ONTOLOGY GSEA (BIOLOGICAL PROCESS, CELLULAR COMPONENT, MOLECULAR FUNCTION)
set.seed(20240306)
df_gsea_res <- map_dfr(
    .x = c("BP", "CC", "MF"),
    .f = ~ gseGO(ranked_gene_list,
                 ont = .x,
                 keyType = "ENSEMBL",
                 OrgDb = "org.Hs.eg.db",
                 pvalueCutoff = 1,
                 eps = 1e-300
    ) %>% 
        as_tibble() %>% 
        clean_names %>% 
        mutate(ontology = .x, .before = 1)
)
df_gsea_res

# plot
df_gsea_res %>% 
    filter(p_adjust < 0.05) %>% 
    mutate(enrichment_sign = ifelse(nes < 0, "neg", "pos")) %>% 
    group_by(ontology, enrichment_sign) %>% 
    mutate(order = row_number(),
           label = ifelse(order %in% seq(1, 5) | -log10(p_adjust) > 12, description %>% str_wrap(30), NA)
    ) %>% 
    
    # plot
    ggplot(aes(x = -log10(p_adjust), y = nes)) +
    geom_point(aes(size = set_size, fill = ontology, alpha = abs(nes)), color = "black", shape = 21) +
    geom_vline(aes(xintercept = -log10(0.05)), color = "black") +
    geom_text_repel(aes(label = label), min.segment.length = 0, size = 3,
                    box.padding = 0.75, max.overlaps = 40) +
    facet_wrap(vars(ontology), nrow = 1) +
    scale_size_continuous(range = c(1, 5), limits = c(10, 500)) +
    scale_alpha_continuous(range = c(0.1, 1), guide = "none") +
    scale_fill_manual(values = ontology_colors, guide = "none") +
    scale_color_manual(values = ontology_colors, guide = "none") +
    guides(size = guide_legend(title = "n genes")) +
    labs(y = "Normalized enrichment score", x = "-log10(FDR)", 
         title = "Gene set ontology enrichments"
    ) +
    theme(legend.position = c(0.97, 0.50),
          legend.key.size = unit(0.2, "cm")
    )

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/5G.GSEA_biological", .x), # or grcca/withGray
                  width = 16, height = 4.5)
)

# export table
write_xlsx(list("Biological process" = as.data.frame(df_gsea_res %>% filter(ontology == "BP") %>% dplyr::select(-ontology)),
                "Cellular component" = as.data.frame(df_gsea_res %>% filter(ontology == "CC") %>% dplyr::select(-ontology)),
                "Molecular function" = as.data.frame(df_gsea_res %>% filter(ontology == "MF") %>% dplyr::select(-ontology))
),
paste0(base_dir, "outputs/tables/for_manuscript/TableS8_TRANSCRIPTS_GSEA_res.xlsx")
)



# DEPR GSEA ---------------------------------------------------------------



######## ******* EXTRAS **********


# F | CELL TYPE  
m_cell <- read_csv(paste0(base_dir, "data/cell_genes_lake++ (1).csv")) %>% 
    left_join(df_hsapiens_genome %>% 
                  dplyr::select(gene_name, gene_id) %>% 
                  dplyr::rename("gene" = "gene_name"),
              by = join_by(gene)
    ) %>% 
    dplyr::select(label, gene_id) %>% 
    dplyr::rename_all(~c("gs_name", "ensembl_gene"))

set.seed(20240306)
df_cell_res <- GSEA(ranked_gene_list, TERM2GENE = m_cell, pvalueCutoff = 1) %>% as_tibble() %>% clean_names()

df_cell_res %>% 
    mutate(label = ifelse(p_adjust < 0.05, description, NA)) %>% 
    
    ggplot(aes(x = -log10(p_adjust), y = nes)) +
    geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
    geom_point(aes(size = set_size, fill = nes), shape = 21) +
    geom_text_repel(aes(label = label), min.segment.length = 0, size = 3,
                    box.padding = 0.75, max.overlaps = 10) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, guide = "none") +
    scale_size_continuous(range = c(3, 6), limits = c(10, 500)) +
    #scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
    guides(size = guide_legend(title = "n genes")) +
    labs(y = "Normalized enrichment score", x = "-log10(FDR)", title = bquote(bold("F")~" | Gene set cell-type enrichments")) +
    theme(legend.position = "right")

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/GSEA_cell_type", .x), # or grcca/withGray
                  width = 4.5, height = 1.5)
)

# G | POSITIONAL
m_position <- msigdbr(species = "Homo sapiens", category = "C1") %>%
    dplyr::select(gs_name, ensembl_gene)

sig_pge_res <- GSEA(sig_gene_list, TERM2GENE = m_position, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
head(sig_pge_res)

# --> position hypergeometric test
gene_universe <- df_x_res %>% pull(ensembl_gene_id) %>% unique
df_sig_chr_count <- m_position %>% 
    filter(ensembl_gene %in% names(sig_gene_list)) %>% 
    left_join(df_x_res %>% 
                  filter(p_adj < 0.05 & abs(z_score) >= 2) %>% 
                  dplyr::select(ensembl_gene_id, ensembl_transcript_id),
              by = join_by("ensembl_gene" == "ensembl_gene_id")
    ) %>% 
    count(gs_name, sort = TRUE)

df_chr_overlap <- tibble()
for (chr in unique(df_sig_chr_count$gs_name))  {
    
    chr_genes <- m_position %>% 
        filter(gs_name == chr) %>% 
        pull(ensembl_gene)
    
    chr_grcca_intersect <- intersect(names(sig_gene_list), chr_genes)
    
    p <- 1 - phyper(q = length(chr_grcca_intersect), 
                    m = length(sig_gene_list), 
                    n = length(gene_universe) - length(sig_gene_list), 
                    k = length(chr_genes)
    )
    
    df_tmp <- tibble(position = chr,
                     overlap = length(chr_grcca_intersect),
                     n_position_genes = length(chr_genes),
                     p_value = p)
    
    df_chr_overlap <- df_chr_overlap %>% bind_rows(df_tmp)
    
}
#plot hypergeometric test
df_chr_overlap %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>% 
    arrange(p_value) %>% 
    filter(p_value < 0.05) %>% 
    
    ggplot(aes(x = reorder(position, p_value), y = -log10(p_value))) +
    geom_point() +
    geom_hline(aes(yintercept = -log10(0.05)))

# plot n
df_sig_chr_count %>% 
    filter(n > 2) %>% 
    ggplot(aes(x = reorder(gs_name, -n), y = n)) +
    geom_col(aes(fill = n), color = "black") +
    scale_fill_gradient(low = "white", high = "maroon", limits = c(0, 8), guide = "none") +
    labs(x = "chromosome position", n = "number of significant transcripts in position",
         title = bquote(bold("G")~" | Chromosome position")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/chromsome_position_count", .x), # or grcca/withGray
                  width = 8, height = 4)
)

# H | BIOTYPE
sig_transcript_list <- df_x_res %>%
    filter(p_adj < 0.05 & abs(z_score) >= 2) %>% 
    #filter(abs(z_score) >= 3) %>% 
    dplyr::select(ensembl_transcript_id, pearsons_r) %>% 
    arrange(-pearsons_r) %>%
    deframe # n = 250

m_biotype <- df_hsapiens_genome %>% 
    dplyr::select(transcript_biotype, transcript_id) %>% 
    dplyr::rename_all(~ c("gs_name", "ensembl_transcript"))

sig_biotype_res <- GSEA(sig_transcript_list, TERM2GENE = m_biotype, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
head(sig_biotype_res)

sig_biotype_res %>% 
    filter(pvalue < 0.05) %>% 
    mutate(description = str_replace_all(description, "_", " ")) %>% 
    plot_gsea_res(plot.title = bquote(bold("H")~" | Biotype"), 
                  y.limits = c(0, 7), size.limits = c(10, 115), fill.scale.high = "maroon") +
    theme(legend.position = "none")

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/GSEA_biotype", .x), # or grcca/withGray
                  width = 4, height = 2)
)


# SAVE TABLES FOR PUBLICATION

write_xlsx(list("GO functional GSEA" = sig_go_res,
                "Position GSEA" = sig_pge_res,
                "Cell type GSEA" = sig_cell_res,
                "Biotype GSEA" = sig_biotype_res
),
"~/Documents/PhD/manuscripts/BiolPsych2024/tables/TableS8.xlsx"
)



# Fig SXX | Correlate X structure correlations with expression values --------------------------------------

df_x_res %>% filter(p_adj < 0.05 & abs(z_score) >= 2) #n = 131 (_NOGRAY); 257 (_WITHGray)

# ILLUSTRATE STRUCTURE CORRELATION AND EXPRESSION IN SCZ

# normatlized & regressed expression of GRCCA Pearson's r
df_vsd_regress_filt %>% 
    pivot_longer(contains("ENST"), names_to = "ensembl_transcript_id", values_to = "expression_value") %>% 
    mutate(dx = case_when(
        str_detect(sample, "control") ~ "Control",
        str_detect(sample, "bipolar") ~ "BD",
        str_detect(sample, "mdd") ~ "MDD",
        str_detect(sample, "schizo") ~ "SCZ"
    ) %>% factor(levels = names(dx_colors))) %>% 
    group_by(ensembl_transcript_id, dx) %>% 
    summarise(mean_expr = mean(expression_value)) %>% 
    left_join(df_x_res) %>% 

    ggplot(aes(x = mean_expr, y = pearsons_r)) +
    geom_hline(yintercept = 0, lty = 2, color = "black") +
    geom_vline(xintercept = 0, lty = 2, color = "black") +
    geom_point(aes(color = pearsons_r, alpha = abs(pearsons_r))) +
    geom_smooth(method = "lm", color = "black") +
    stat_cor(aes(label = after_stat(r.label)), label.sep = "\n", label.y.npc = "bottom", vjust = -0.8) +
    facet_wrap(vars(dx), nrow = 1) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    labs(x = "Mean normalized & corrected expression", y = "Structure correlation (r)",
         title = bquote(bold("C")~" | Correlations of GRCCA structure correlations with mean expression value across all transcripts")) +
    theme(legend.position = "none")

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/xres_struc_cor_mean_expr", .x), # or grcca/withGray
                  width = 12, height = 4)
)







# Fig SXX | compare with DETs -------------------------------------------------------

# ANALYSIS

# load results from Nirmala's DGE analysis
load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res

# identify DETs
DETs <- df_akula_res %>% 
  filter(level == "transcript" & padj < 0.05) %>% 
  pull(id) %>% 
  unique # 1008

# identify significant GRCCA transcripts
grcca_transcripts <- df_x_res %>%
  filter(p_adj < 0.05) %>%
  pull(ensembl_transcript_id) # 120

# identify GRCCA hits that are also differentially expressed
grcca_DETs <- intersect(DETs, grcca_transcripts) # 1

df_modules_filt %>% filter(ensembl_transcript_id == "ENST00000376762")


# total number of transcripts tested in GRCCA
total_transcripts <- df_x_res %>% pull(ensembl_transcript_id) %>% unique

# hypergeometric overrepresentation of DETs by module

# q = number of DETs in the module
# m = total DETs in universe
# n = transcripts in universe that were *not* DETs (total transcripts - total DETs)
# k = total transcripts in module

df_hypergeometric_dets <- df_modules_filt %>% 
  filter(ensembl_transcript_id %in% DETs) %>% 
  dplyr::count(module, color) %>% 
  dplyr::rename("q" = "n") %>% 
  left_join(
    df_modules_filt %>% 
      dplyr::count(module, color) %>% 
      dplyr::rename("k" = "n")
  ) %>% 
  mutate(m = length(DETs),
         n = length(total_transcripts) - length(DETs),
         p_value = phyper(q, m, n, k, lower.tail = FALSE),
         p_adj = p.adjust(p_value, method = "fdr")
  ) %>% 
  arrange(p_value)

# PLOT
df_hypergeometric_dets %>%
  ggplot() +
  geom_segment(aes(x = module, xend = module, y = 0, yend = -log10(p_adj))) +
  #geom_col(width = 0.05, fill = "black") +
  geom_point(aes(x = module, y = -log10(p_adj), fill = I(color), size = (q/k)*100), shape = 21) +
  geom_hline(aes(yintercept = -log10(0.05)), lty = 2) +
  guides(size = guide_legend(title = "Percent of module that was \n differentially expressed")) +
  labs(x = "Module", y = "-log10(hypergeometric FDR)",
       title = "Modules with overrepresentation of DETs") +
  theme(legend.position = c(0.15, 0.80))

# save to project dir
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/DETs", .x), # or grcca/withGray
                width = 7, height = 4)
)


# save for publication
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0("~/Documents/PhD/manuscripts/BiolPsych2024/figures/supplement/FigS7", .x),
                width = 7, height = 6)
)

# save table
df_tableS7 <- df_akula_res %>% 
  filter(id %in% grcca_transcripts & padj < 0.05) %>% 
  left_join(df_modules_filt %>% dplyr::rename("id" = "ensembl_transcript_id")) %>% 
  dplyr::select(-c(level, ensembl_gene_id, gene_symbol)) %>% 
  arrange(comparison, module, padj) %>% 
  as.data.frame
write.xlsx(df_tableS7, file = "~/Documents/PhD/manuscripts/BiolPsych2024/tables/TableS7.xlsx", sheetName = "DETs", row.names = FALSE)


