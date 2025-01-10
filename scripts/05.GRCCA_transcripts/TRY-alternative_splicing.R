
################################################################################

# Characterize role of alternative splicing in GRCCA results
# Distinguish results from gene- vs transcript-level

################################################################################


base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/full_analysis_scripts/transcripts/overall_subset/GRCCA/04.3.load_GRCCA_results.R"))

# Load gene & transcript-level GRCCA results --------------------------------------------------------------------

# LOAD GENE GRCCA RESULTS (X RES ONLY)
df_x_res_genes <- read_excel(paste0(base_dir, "outputs/tables/for_manuscript/TableS3_GENES_GRCCA_res_withGray.xlsx"), sheet = 3) %>% 
    mutate(module = as.numeric(module %>% str_remove("geneM")),
           module = factor(paste0("geneM", module), levels = paste0("geneM", 0:max(module)))
           )

# LOAD TRANSCRIPT GRCCA RESULTS
df_y_res_transcripts <- read_excel(paste0(base_dir, "outputs/tables/for_manuscript/TableS7_TRANSCRIPTS_GRCCA_res_withGray.xlsx"), sheet = 2)
df_x_res_transcripts <- read_xlsx(paste0(base_dir, "outputs/tables/for_manuscript/TableS7_TRANSCRIPTS_GRCCA_res_withGray.xlsx"), sheet = 3) %>% 
    mutate(module = as.numeric(module %>% str_remove("transcriptM")),
           module = factor(paste0("transcriptM", module), levels = paste0("transcriptM", 0:max(module)))
    )

# TRANSCRIPT RAW COUNTS
df_transcript_raw_counts <- read.csv(paste0(base_dir, "data/transcript_count_matrix.csv")) %>% 
  as_tibble() %>% 
  dplyr::rename("ensembl_transcript_id" = "X") %>% 
  rename_all(~str_remove(.x, "_MergedBam_stringtieOutput")) %>% 
  clean_names() %>% 
  rename_all(~str_remove(.x, "x")) 


# Identify rare transcripts -----------------------------------------------

# DEFINE TRANSCRIPT UNIVERSE AS ALL TRANSCRIPTS IN GRCCA ANALYSIS
transcript_universe <- df_x_res_transcripts %>% pull(ensembl_transcript_id) %>% unique
length(transcript_universe) # n = 54302

# identify rare transcripts (at least 10 raw counts in at least 80% of samples)
rare_transcripts <- df_transcript_raw_counts %>% 
  filter(ensembl_transcript_id %in% transcript_universe) %>% 
  mutate_if(is.numeric, ~ifelse(.x >= 10 & .x < 100, 1, 0)) %>% 
  mutate(thresh = rowSums(dplyr::select(., -ensembl_transcript_id))/(ncol(.) - 1), .before = 2) %>% 
  filter(thresh >= 0.80) %>% 
  pull(ensembl_transcript_id) 
length(rare_transcripts) # n = 23211

# add rare column to x res for analysis
df_x_res_transcripts <- df_x_res_transcripts %>% 
  mutate(prevalence = ifelse(ensembl_transcript_id %in% rare_transcripts, "rare", "common"))


# Identify genes represented in isoform-level analysis, and number of transcripts associated with each ---------------------

# KNOWN GENES REPRESENTED BY TRANSCRIPTS IN GRCCA ANALYSIS
gene_universe <- df_x_res_transcripts %>%
    filter(!is.na(ensembl_gene_id)) %>% 
    pull(ensembl_gene_id) %>% 
    unique()
length(gene_universe) # n = 18031

# COUNT NUMBER OF TRANSCRIPTS PER GENE IN ISOFORM-LEVEL GRCCA ANALYSIS
df_number_of_transcripts_per_gene <- df_x_res_transcripts %>% 
    dplyr::select(ensembl_gene_id, gene_symbol, ensembl_transcript_id) %>% 
    distinct() %>% 
    dplyr::count(ensembl_gene_id, gene_symbol, sort = TRUE) %>% 
    filter(!is.na(ensembl_gene_id)) # 18031 unique genes

# GENES WITH MULTIPLE TRANSCRIPTS IN GRCCA ANALYSIS
genes_with_multiple_transcripts <- df_number_of_transcripts_per_gene %>% filter(n > 1) %>% pull(ensembl_gene_id) %>% unique
length(genes_with_multiple_transcripts) # 9996

# GENES WITH ONLY ONE TRANSCRIPT IN GRCCA ANALYSIS
genes_with_one_transcript <- df_number_of_transcripts_per_gene %>% filter(n == 1) %>% pull(ensembl_gene_id) %>% unique
length(genes_with_one_transcript) # 8035


# Identify genes that were implicated in isoform-level analysis, but not gene-level ---------------------

# GENE-LEVEL GRCCA SIGNIFICANT GENES
grcca_genes_gene_level <- df_x_res_genes %>% 
    filter(significant == 1) %>% 
    pull(ensembl_gene_id) %>% 
    unique()
length(grcca_genes_gene_level) # 2026

# ISOFORM-LEVEL GRCCA SIGNIFICANT GENES (at least 1 isoform significant)
grcca_transcripts <- df_x_res_transcripts %>% 
    filter(significant == 1) %>% 
    pull(ensembl_transcript_id) %>% 
    unique
df_sig_transcripts_per_gene <- df_x_res_transcripts %>% 
    filter(ensembl_transcript_id %in% grcca_transcripts & !is.na(ensembl_gene_id)) %>% 
    count(ensembl_gene_id, gene_symbol, sort = TRUE)

# pull IDs
grcca_genes_transcript_level <- df_sig_transcripts_per_gene %>% 
    pull(ensembl_gene_id) %>% 
    unique()
length(grcca_genes_transcript_level) # 234
length(grcca_transcripts) # 257



# Fig 6A | Number of genes implicated by transcript-level GRCCA that were not identified at gene level --------------------

# LIST GENES THAT WERE IMPLICATED BOTH AT GENE AND TRANSCRIPT LEVEL 
grcca_genes_both_levels <- intersect(grcca_genes_gene_level, grcca_genes_transcript_level)
length(grcca_genes_both_levels) # 100

# Genes identified at transcript level but *not* gene level
df_x_res_transcripts_unique <- df_x_res_transcripts %>% 
    filter(significant == 1) %>% 
    filter(!(ensembl_gene_id %in% grcca_genes_both_levels)) # 134 unique known genes identified by transcript analysis (map to 146 transcripts)

## PLOT (Venn diagram)


## WRITE TABLE
df_x_res_transcripts_unique %>% 
    write_xlsx(path = paste0(base_dir, "outputs/tables/for_manuscript/TableS9_transcript_level_only_results.xlsx"))
  

# Fig S8 | Genes implicated both at gene and transcript level ------------


# PLOT
df_x_res_transcripts %>% 
    filter((ensembl_gene_id %in% grcca_genes_both_levels) & significant == 1) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), transcript_symbol, gene_symbol)) %>% 
    
    ggplot(aes(x = pearsons_r, y = reorder(gene_symbol, pearsons_r))) +
    geom_point(aes(color = module)) +
    scale_color_manual(values = module_colors) +
    labs(x = "Pearson's r", y = "",
         title = "Genes implicated by gene- and transcript-level analyses")



# Fig S9 | Transcripts not implicated at gene level (not included for now)----------------------------------------------------------------


# PLOT: TRANSCRIPT-LEVEL SPECIFIC GENES
df_x_res_transcripts_unique %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), transcript_symbol, gene_symbol)) %>% #count(module, sort = TRUE)
    
    ggplot(aes(x = pearsons_r, y = reorder(gene_symbol, pearsons_r))) +
    geom_point(aes(color = module)) +
    scale_color_manual(values = module_colors) +
    labs(x = "Pearson's r", y = "",
         title = "Genes implicated by transcript-level analysis only")



# Fig SXX | Heatmap of transcript-level only transcripts ------------------


## GET EXPERSSION VALUES FOR UNIQUE TRANSCRIPTS & NORMALIZE
df_x_res_transcripts_unique_expr <- df_vsd_regress_filt %>% 
    dplyr::select(sample, all_of(df_x_res_transcripts_unique$ensembl_transcript_id)) %>% 
    # pivot_longer(contains("ENST"), names_to = "ensembl_transcript_id", values_to = "expr") %>% 
    # # noramlize within transcript
    # group_by(ensembl_transcript_id) %>% 
    # mutate(norm_expr = scale(expr)[,1]) %>% 
    
    mutate_if(is.numeric, ~ scale(.x)[,1]) %>% 
    pivot_longer(contains("ENST"), names_to = "ensembl_transcript_id", values_to = "norm_expr") %>% 
    
    # add transcript symbol
    left_join(df_ensembl_to_symbol) %>% 
    mutate(transcript_symbol = ifelse(is.na(transcript_symbol), ensembl_transcript_id, transcript_symbol))
    
## CREATE MATRIX FOR PLOTTING ORDER
m_unique_transcripts_expr <- df_x_res_transcripts_unique_expr %>% 
    pivot_wider(id_cols = sample, names_from = transcript_symbol, values_from = norm_expr) %>% 
    column_to_rownames("sample") %>% 
    as.matrix
sample_order <- m_unique_transcripts_expr[hclust(dist(m_unique_transcripts_expr))$order,] %>% rownames
transcript_order <- m_unique_transcripts_expr[hclust(dist(m_unique_transcripts_expr))$order,] %>% colnames

transcript_order <- df_x_res_transcripts_unique %>% 
    arrange(pearsons_r) %>% 
    pull(transcript_symbol)

pheatmap(m_unique_transcripts_expr)

## PLOT
df_x_res_transcripts_unique_expr %>% 
    mutate(dx = str_remove(sample, ".*_"),
           dx = case_when(
               dx == "control" ~ "Control",
               dx == "bipolar" ~ "BD",
               dx == "mdd" ~ "MDD",
               dx == "schizo" ~ "SCZ"
           ) %>% factor(levels = names(dx_colors))
    ) %>% 
    arrange(dx) %>% 
    mutate(sample = factor(sample, levels = sample_order)) %>% 
    mutate(transcript_symbol = factor(transcript_symbol, levels = transcript_order)) %>% 

    ggplot(aes(x = sample, y = transcript_symbol)) +
    geom_tile(aes(fill = norm_expr)) +
    scale_fill_viridis() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))




# Fig 6B | Example expression patterns of transcripts identified in transcript-level analysis only -----------------------------------------


# IDENTIFY GENES ASSOCIATED WITH TRANSCRIPTS THAT WERE UNIQUELY IDENTIFIED IN TRANSCRIPT-LEVEL ANALYSIS
df_x_res_transcripts %>% 
    filter(ensembl_transcript_id %in% tx_level_only_ids) %>% filter(is.na(ensembl_gene_id)) # 7 have unknown gene ID
tx_level_only_GENEids <- df_x_res_transcripts %>% 
    filter(ensembl_transcript_id %in% tx_level_only_ids & !is.na(ensembl_gene_id)) %>% # 139 have known gene ID
    pull(ensembl_gene_id) %>% 
    unique
length(tx_level_only_GENEids) # n = 134


# PLOT    
df_x_res_transcripts %>% 
    
    # select genes that had txt-level only significant transcripts
    filter(ensembl_gene_id %in% tx_level_only_GENEids) %>% 
    
    # label significant transcripts with crosses?
    mutate(transcript_symbol =  str_remove(transcript_symbol, ".*-"),
           label = ifelse(significant == 1, "+", "")
    ) %>% 
    
    # add expression data
    left_join(df_vsd_regress_filt %>% 
                  pivot_longer(contains("ENST"), names_to = "ensembl_transcript_id", values_to = "expr"),
              by = join_by(ensembl_transcript_id)
    ) %>% 
    
    # add diagnosis information
    left_join(df_covariates %>% dplyr::select(sample, dx)) %>% 
    mutate(dx = factor(dx, levels = names(dx_colors))) %>% 
    
    # filter for just SCZ vs controls
    filter(dx %in% c("SCZ", "Control")) %>% 
    filter(!is.na(sample)) %>% 
    
    # SELECT GENE HERE
    #filter(gene_symbol == "WDFY4") %>% 
    #filter(gene_symbol == "KIAA2026") %>% 
    #filter(gene_symbol == "CHGA") %>% 
    filter(gene_symbol %in% c("ARFGEF1", "CHGA", "WDFY4")) %>% 
    
    # plot!
    ggplot(aes(x = transcript_symbol, y = expr)) +
    geom_point(aes(fill = dx), shape = 21, size = 1.25, 
               position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75)) +
    geom_boxplot(aes(color = dx), fill = "transparent", outlier.shape = NA) +
    geom_text(aes(label = label, y = -Inf), check_overlap = TRUE, vjust = -2, size = 7) +
    facet_wrap(vars(gene_symbol), scales = "free") +
    force_panelsizes(cols = c(4, 5, 2)) +
    scale_fill_manual(values = dx_colors) +
    scale_color_manual(values = dx_colors) +
    guides(fill = guide_legend(title = ""), color = guide_legend(title = "")) +
    labs(x = "Isoform number", y = "Normalized & corrected expression value",
         caption = "+ significant isoform",
         title = "Isoform expression in SCZ vs control"
    ) +
    theme(legend.position = "bottom")

#df_x_res_genes %>% filter(gene_symbol == "ARFGEF1")

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/6B.tx_level_only_transcripts_examples", .x), # or grcca/withGray
                  width = 8, height = 4)
)



# Fig 6C | Plot retained intron example ---------------------------------------


## IDENTIFY BIOTYPE FOR EACH TRANSCRIPT ID
df_transcript_biotype <- df_hsapiens_genome %>% 
    filter(type == "transcript") %>% 
    dplyr::select(transcript_id, transcript_biotype)

# IDENTIFY RETAINED INTRONS IN RESULTS
df_transcript_biotype %>% 
    filter(transcript_biotype == "retained_intron" & transcript_id %in% tx_level_only_ids) %>% 
    left_join(df_ensembl_to_symbol, by = join_by("transcript_id" == "ensembl_transcript_id")) %>% arrange(transcript_symbol) %>% print(n = nrow(.))

# IDENTIFY GENES ASSOCIATED WITH TRANSCRIPTS THAT WERE UNIQUELY IDENTIFIED IN TRANSCRIPT-LEVEL ANALYSIS
tx_level_only_GENEsymbols <- df_x_res_transcripts %>% 
    filter(ensembl_transcript_id %in% tx_level_only_ids & !is.na(ensembl_gene_id)) %>% # 139 have known gene ID
    arrange(-abs(pearsons_r)) %>% 
    pull(gene_symbol) %>% 
    unique
length(tx_level_only_GENEsymbols) # n = 134

# IDENTIFY EXONS FOR SIGNIFICANT TRANSCRIPTS
exons <- df_hsapiens_genome %>% 
    filter((gene_name %in% tx_level_only_GENEsymbols) & (transcript_id %in% transcript_universe)) %>% 
    mutate(transcript_name = ifelse(is.na(transcript_name), transcript_id, transcript_name)) %>% 
    filter(!is.na(transcript_name)) %>% 
    mutate(significant = ifelse(transcript_id %in% tx_level_only_ids, "+", "")) %>% 
    filter(type == "exon") %>% 
    mutate(start = start/1e6,
           end = end/1e6) %>% 
    mutate(transcript_name = str_remove(transcript_name, ".*-"),
           transcript_biotype = str_replace_all(transcript_biotype, "_", " "),
           gene_name = factor(gene_name, levels = tx_level_only_GENEsymbols)
    )

# example plot (main text)
exons %>% 
    
    # select gene to plot isoforms of
    filter(gene_name == "CHGA") %>% 
    
    # plot isoforms
    ggplot(aes(
        xstart = start,
        xend = end,
        y = transcript_name
    )) +
    geom_range(
        aes(fill = transcript_biotype)
    ) +
    geom_intron(
        data = to_intron(exons %>% filter(gene_name == "CHGA"), "transcript_name"),
        aes(strand = strand)
    ) +
    geom_text(aes(label = significant, x = -Inf), hjust = -1, size = 7, check_overlap = "TRUE") +
    facet_wrap(vars(gene_name)) +
    guides(fill = guide_legend(title = "Biotype", title.position = "top", title.hjust = 0.5)) +
    labs(y = "Isoform number", x = "Position (Mb)", caption = "+ significant isoform",
         title = "CHGA isoform biotypes"
    ) +
    theme(legend.position = "bottom")

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/6C.tx_level_only_transcripts_CHGA", .x), # or grcca/withGray
                  width = 8, height = 4)
)






# Fig 6D | Transcripts not implicated at gene level: hypergeometric module test ---------------------------
    

## HYPERGEOMETRIC TESTING      
    
    # hypergeometric overrepresentation of tx-level transcripts by module
    # pool = transcripts in analysis
    
    # pull IDs for trancsript-level only IDs
    tx_level_only_ids <- df_x_res_transcripts_unique %>% 
        #filter(ensembl_gene_id %in% df_x_res_genes$ensembl_gene_id) %>% # turn on to include transcripts from genes included in gene-level analysis
        pull(ensembl_transcript_id) %>% 
        unique

    # run test
    df_module_overlap <- tibble()
    for (mod in unique(df_x_res_transcripts$module))  {
        
        mod_transcripts <- df_x_res_transcripts %>% 
            filter(module == mod) %>% 
            pull(ensembl_transcript_id)
        
        mod_grcca_intersect <- intersect(tx_level_only_ids, mod_transcripts)
        
        p <- 1 - phyper(q = length(mod_grcca_intersect), 
                        m = length(tx_level_only_ids), 
                        n = length(transcript_universe) - length(tx_level_only_ids), 
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
        arrange(p_adj) %>% 
        mutate(module = factor(module, levels = levels(df_modules_filt$module))) %>% 
        mutate(p_adj = ifelse(overlap == 0, 1, p_adj))

# PLOT RES
    df_module_overlap %>%
        mutate(label = ifelse(p_adj < 0.01, paste0("n = ", overlap), NA)) %>% 
        
        ggplot() +
        geom_segment(aes(x = module, xend = module, y = 0, yend = -log10(p_adj))) +
        #geom_col(width = 0.05, fill = "black") +
        geom_point(aes(x = module, y = -log10(p_adj), fill = module, size = (overlap/mod_size)*100), shape = 21) +
        geom_text(aes(x = module, y = -log10(p_adj), label = label), vjust = 0.5, hjust = -0.25) +
        geom_hline(aes(yintercept = -log10(0.01)), linewidth = 0.25) +
        annotate(geom = "text", label = "FDR < 0.01", x = 1, y = -log10(0.01), hjust = 0) +
        scale_fill_manual(values = module_colors, guide = "none") +
        guides(size = guide_legend(title = "Percent of module that was \n significant")) +
        ylim(c(0, 16)) +
        coord_flip() +
        labs(x = "", y = "-log10(hypergeometric FDR)",
             title ="Module overrepresentation of \ntranscript-level only results") +
        theme(legend.position = c(0.75, 0.85)
              )
    
    # save
    map(
        .x = c(".png", ".pdf"),
        .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/6D.tx_level_only_transcripts_hypergeometric", .x), # or grcca/withGray
                      width = 5, height = 6)
    )

    
# Fig 6E | cell type hypergeometric ---------------------------------------------
    
    ## LOAD CELL TYPE DATA
    load(paste0(base_dir, "objects/cell_type_data.RDS")) # df_cell_type
    
    ## RUN HYPERGEOMETRIC TEST ON TX-LEVEL ONLY TRANSCRIPTS FOR EACH CELL TYPE  
    gene_universe <- df_x_res_transcripts %>% 
        pull(ensembl_gene_id) %>% 
        unique
    cell_types <- df_cell_type %>% filter(str_detect(type, "lake")) %>% pull(type) %>% unique
    
    df_cell_type_enrichment <- tibble()
    for (i in cell_types) {
        
        cell_type_genes <- df_cell_type %>% 
            filter(type == i & ensembl_gene_id %in% gene_universe & !(is.na(ensembl_gene_id))) %>% 
            pull(ensembl_gene_id)
        cell_type_genes_sig_transcripts <- intersect(tx_level_only_GENEids, cell_type_genes)
        
        p <- 1 - phyper(q = length(cell_type_genes_sig_transcripts), # number of significant genes that are of that cell type
                        m = length(cell_type_genes), # number of genes that are of that cell type total
                        n = length(gene_universe) - length(cell_type_genes), # number of genes that are not of that cell type total
                        k = length(tx_level_only_GENEids) # number of significant genes
        )
        
        df_tmp <- tibble(
            cell_type = i,
            overlap = length(cell_type_genes_sig_transcripts),
            phyper = p
        )
        df_cell_type_enrichment <- df_cell_type_enrichment %>% bind_rows(df_tmp)
        
    }
    
    df_cell_type_enrichment <- df_cell_type_enrichment %>% 
        mutate(p_adj = p.adjust(phyper, method = "fdr")) %>% 
        arrange(p_adj)
    df_cell_type_enrichment
    
    # PLOT
    df_cell_type_enrichment %>% 
        mutate(cell_type = case_when(
            cell_type == "neuro-ex_lake" ~ "Neuro-Ex",
            cell_type == "neuro-in_lake" ~ "Neuro-In",
            cell_type == "endo_lake" ~ "Endo",
            cell_type == "synapses_lake" ~ "Synapses",
            cell_type == "micro_lake" ~ "Micro",
            cell_type == "opc_lake" ~ "OPC",
            cell_type == "oligo_lake" ~ "Oligo",
            cell_type == "astro_lake" ~ "Astro",
        )
        ) %>% 
        mutate(label = paste0("n = ", overlap)) %>% 
        
        ggplot(aes(x = reorder(cell_type, p_adj), y = -log10(p_adj))) +
        geom_col(aes(fill = -log10(p_adj))) +
        geom_hline(aes(yintercept = -log10(0.05)), lty = 2) +
        geom_text(aes(label = label), vjust = -0.25) +
        scale_fill_gradient(low = "white", high = "#222222", limits = c(0, 3.2)) +
        labs(x = "", y = "-log10(hypergeometric FDR)",
             title = "Overrepresented cell types"
        ) +
        theme_classic() +
        theme(axis.text = element_text(size = 10),
              legend.position = "none")
    
    # save
    map(
        .x = c(".png", ".pdf"),
        .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/6E.tx_level_only_transcripts_cell_type", .x), # or grcca/withGray
                      width = 5.5, height = 4.5)
    )
    
    
# Fig 6F | Biotype enrichment -------------------------------------------------
    
    
    ## RUN HYPERGEOMETRIC TEST ON TX-LEVEL ONLY TRANSCRIPTS FOR EACH BIOTYPE  
    biotypes <- df_transcript_biotype %>% pull(transcript_biotype) %>% unique
    
    df_biotype_enrichment <- tibble()
    for (i in biotypes) {
        
        biotype_transcripts <- df_transcript_biotype %>% 
            filter(transcript_biotype == i & transcript_id %in% transcript_universe & !(is.na(transcript_id))) %>% 
            pull(transcript_id)
        biotype_transcripts_sig_transcripts <- intersect(tx_level_only_ids, biotype_transcripts)
        
        p <- 1 - phyper(q = length(biotype_transcripts_sig_transcripts), # number of significant genes that are of that cell type
                        m = length(biotype_transcripts), # number of genes that are of that cell type total
                        n = length(transcript_universe) - length(cell_type_genes), # number of genes that are not of that cell type total
                        k = length(tx_level_only_ids) # number of significant genes
        )
        
        df_tmp <- tibble(
            biotype = i,
            overlap = length(biotype_transcripts_sig_transcripts),
            phyper = p
        )
        df_biotype_enrichment <- df_biotype_enrichment %>% bind_rows(df_tmp)
        
    }
    
    df_biotype_enrichment <- df_biotype_enrichment %>% 
        mutate(p_adj = p.adjust(phyper, method = "fdr")) %>% 
        arrange(p_adj)
    df_biotype_enrichment %>% filter(overlap != 0)
    
    # PLOT
    df_biotype_enrichment %>% 
        filter(overlap != 0) %>% 
        mutate(biotype = str_replace_all(biotype, "_", " ") %>% str_wrap(width = 10)) %>% 
        mutate(label = paste0("n = ", overlap)) %>% 
        
        ggplot(aes(x = reorder(biotype, p_adj), y = -log10(p_adj))) +
        geom_col(aes(fill = -log10(p_adj))) +
        geom_hline(aes(yintercept = -log10(0.05)), lty = 2) +
        geom_text(aes(label = label), vjust = -0.25) +
        scale_fill_gradient(low = "white", high = "#222222", limits = c(0, 3)) +
        labs(x = "", y = "-log10(hypergeometric FDR)",
             title = "Overrepresented biotypes"
        ) +
        theme_classic() +
        theme(axis.text = element_text(size = 10),
              legend.position = "none")
    
    # save
    map(
        .x = c(".png", ".pdf"),
        .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/6F.tx_level_only_transcripts_biotype", .x), # or grcca/withGray
                      width = 5.5, height = 4.5)
    )
    
    
    

# Fig 6G | GO enrichment of tx-level only transcripts -------------------------


# map each gene (that maps to multiple transcripts) to GO term
    go_res <- enrichGO(gene = tx_level_only_GENEids,
                       OrgDb = "org.Hs.eg.db",
                       universe = df_x_res_transcripts %>% pull(ensembl_gene_id) %>% unique,
                       keyType = "ENSEMBL",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       readable = TRUE
    ) %>% as_tibble() %>% clean_names
    
    
# plot
    go_res %>% 
        mutate(`survives FDR` = ifelse(p_adjust < 0.05, "yes", "no")) %>% 
        mutate(description = str_wrap(description, width = 50)) %>%
        filter(pvalue < 0.05) %>% 
        group_by(ontology) %>% top_n(10, wt = -pvalue) %>% 
        mutate(gene_ratio = eval(parse(text = gene_ratio))) %>% 
        #filter(p_adjust < 0.05) %>% 
        
        # plot
        ggplot(aes(x = -log10(pvalue), y = reorder(description, -pvalue))) +
        #y = reorder_within(str_wrap(description, width = 35), -pvalue, module))) +
        geom_point(aes(size = gene_ratio, fill = ontology, color = `survives FDR`),
                   shape = 21, stroke = 1) +
        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
        scale_color_manual(values = c("yes" = "black", "no" = "transparent")) +
        scale_fill_manual(values = ontology_colors, guide = "none") +
        scale_size_continuous(range = c(2, 4)) +
        scale_y_reordered() +
        facet_wrap(vars(ontology), scales = "free_y") +
        labs(y = "", x = "log10(p-value)",
             title = "GO enrichment of transcript-level only results"
        ) +
        guides(size = guide_legend(title = "Gene ratio")) +
        theme(legend.position = "bottom",
              legend.direction = "horizontal",
              legend.box = "horizontal",
              legend.justification = "center",
              axis.text.x = element_text(size = 8),
              legend.key.size = unit(0.3, "cm")
        )
    
    
    map(
        .x = c(".png", ".pdf"),
        .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/6G.tx_level_only_transcripts_GO", .x), # or grcca/withGray
                      width = 13, height = 4)
    )
    
# export results table
    # export table
    write_xlsx(list("Biological process" = as.data.frame(go_res %>% filter(ontology == "BP") %>% dplyr::select(-ontology)),
                    "Cellular component" = as.data.frame(go_res %>% filter(ontology == "CC") %>% dplyr::select(-ontology)),
                    "Molecular function" = as.data.frame(go_res %>% filter(ontology == "MF") %>% dplyr::select(-ontology))
    ),
    paste0(base_dir, "outputs/tables/for_manuscript/TableS10_transcript_level_only_GO_res.xlsx")
    )
    
   
    

# Fig S10 | FOR DISCUSSION: compare immune genes from gene vs transcript level ------------------------

## PULL GSEA RESULTS
    df_gsea_genes <- read_xlsx(paste0(base_dir, "outputs/tables/for_manuscript/TableS4_GENES_GSEA_res.xlsx"), sheet = 1)
    df_gsea_transcripts <- read_xlsx(paste0(base_dir, "outputs/tables/for_manuscript/TableS8_TRANSCRIPTS_GSEA_res.xlsx"), sheet = 1)
    
## IDENTIFY IMMUNE-RELATED GENES IN GENE AND TRANSCRIPT-LEVEL RESULTS
    df_x_res_genes_immune <- df_x_res_genes %>% 
        filter(significant == 1 & module == "geneM6") # geneM6 is the immune module
    df_x_res_transcripts_immune <- df_x_res_transcripts %>% 
        filter(significant == 1 & module %in% c("transcriptM9", "transcriptM18"))
    
## HOW MANY GENES OVERLAP? DO THEY GO THE SAME DIRECTION?
    immune_gene_overlap <- intersect(df_x_res_genes_immune %>% pull(ensembl_gene_id),
              df_x_res_transcripts_immune %>% pull(ensembl_gene_id)
    )
    df_x_res_genes_immune %>% filter(ensembl_gene_id %in% immune_gene_overlap) %>% 
        print(n = nrow(.))
    df_x_res_transcripts_immune %>% filter(ensembl_gene_id %in% immune_gene_overlap) %>% 
        print(n = nrow(.))
    
## PLOT DIRECTION OF IMMUNE-RELATED GENES IMPLICATED AT BOTH GENE AND TRANSCRIPT-LEVEL    
    level_colors <- c("gene" = "#d07444", "transcript1" = "#4C768A", "transcript2" = "#4C768A", "transcript3" = "#4C768A")
    gene_order <- df_x_res_genes_immune %>% filter(ensembl_gene_id %in% immune_gene_overlap) %>% pull(gene_symbol)
    df_x_res_genes_immune %>%
        mutate(level = "gene") %>% 
        bind_rows(
            df_x_res_transcripts_immune %>% 
                group_by(gene_symbol) %>% 
                mutate(transcript_no = row_number()) %>% 
                mutate(level = paste0("transcript", transcript_no))
        ) %>% 
        filter(ensembl_gene_id %in% immune_gene_overlap) %>% 
        mutate(gene_symbol = factor(gene_symbol, levels = gene_order)) %>% 
        
        ggplot(aes(x = pearsons_r, y = gene_symbol)) +
        geom_col(aes(fill = level, alpha = abs(pearsons_r)), color = "black", position = "dodge") +
        scale_fill_manual(values = level_colors) +
        scale_alpha_continuous(limits = c(0, 0.52), guide = "none") +
        labs(x = "Structure correlation", y = "",
             title = "Significant immune gene overlap")
    
    map(
        .x = c(".png", ".pdf"),
        .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/transcripts/grcca/7MCA_in_Ymat_Dim2_in_Cmat/withGray/S10A.immune_gene_intersect_struct_cor", .x), # or grcca/withGray
                      width = 7, height = 7)
    )
    
## PLOT DIRECTION OF IMMUNE-RELATED GENES IMPLICATED AT ONLY GENE- OR TRANSCRIPT-LEVEL
    level_colors <- c("gene" = "#d07444", "transcript" = "#4C768A")
    df_x_res_genes_immune %>%
        mutate(number = 1, level = "gene") %>% 
        filter(!(ensembl_gene_id %in% immune_gene_overlap)) %>% 
        bind_rows(
            df_x_res_transcripts_immune %>% 
                group_by(gene_symbol) %>% 
                mutate(number = row_number()) %>% 
                mutate(level = "transcript") %>% 
                filter(!(ensembl_gene_id %in% immune_gene_overlap))
        ) %>% 

        ggplot(aes(x = pearsons_r, y = reorder_within(gene_symbol, pearsons_r, level))) +
        geom_col(aes(fill = level, alpha = abs(pearsons_r)), color = "black", position = "dodge") +
        facet_wrap(vars(level), scales = "free") +
        scale_fill_manual(values = level_colors) +
        scale_alpha_continuous(limits = c(0, 0.52), guide = "none") +
        scale_y_reordered() +
        labs(x = "Structure correlation", y = "",
             title = "Significant immune gene overlap")
    
    
    
# ### ************************** DEPR ************************* -----------

    
# create table of information ---------------------------------------------


df_x_res %>% 
    filter(transcript_symbol %in% sig_mods_unique_sig_transcripts) %>% 
    arrange(module, pearsons_r) %>% 
    dplyr::select(module, transcript_symbol, pearsons_r, p_adj, prevalence) #%>% 
    
    #write_xlsx(path = paste0(base_dir, "outputs/tables/for_manuscript/unique_sig_mod_transcripts_annotations.xlsx"))

# BIOTYPE
df_hsapiens_genome %>% 
    filter(transcript_name %in% sig_mods_unique_sig_transcripts & type == "transcript") %>% 
    dplyr::select(transcript_name, transcript_biotype)
df_hsapiens_genome %>% filter(transcript_id == "ENST00000435315" & type == "transcript")

# CELL TYPE
load(paste0(base_dir, "objects/cell_type_data.RDS"))

df_cell_type %>% 
    filter(ensembl_gene_id %in% sig_mod_unique_sig_transcripts_geneIDs) %>% 
    print(n = nrow(.))

# POSITION
m_position %>% 
    filter(ensembl_gene %in% sig_mod_unique_sig_transcripts_geneIDs) %>% 
    left_join(df_gene_to_transcript, by = join_by("ensembl_gene" == "ensembl_gene_id")) %>% 
    dplyr::select(gs_name, contains("gene")) %>% 
    distinct()



# literature --------------------------------------------------------------

# Rare coding variants in SCZ
df_x_res %>% filter(gene_symbol == "SETD1A")
df_x_res %>% filter(gene_symbol == "CUL1") # not present in analysis
df_x_res %>% filter(gene_symbol == "XPO7")
df_x_res %>% filter(gene_symbol == "TRIO")
df_x_res %>% filter(gene_symbol == "CACNA1G") # ns
df_x_res %>% filter(gene_symbol == "SP4")
df_x_res %>% filter(gene_symbol == "GRIA3") # ns
df_x_res %>% filter(gene_symbol == "GRIN2A")
df_x_res %>% filter(gene_symbol == "HERC1")
df_x_res %>% filter(gene_symbol == "RB1CC1") # ns

df_x_res %>% filter(gene_symbol == "SMARCE2")

