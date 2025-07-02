
################################################################################

# Response to R3 comment: Transcript-level annotations to increase the novelty
# and interpretability of these findings

################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/Fig5_GRCCAtranscriptsRes/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WTCNA_res.R"))


## Load GRCCA data
load(paste0(base_dir, "objects/GRCCA_transcript_results.Rdata")) # df_lvs, df_y_res, df_x_res



# Load data to annotate individual transcripts ----------------------------



## Look up annotations for transcripts to decide which ones to show
fig5c_sig_genes <- c(
    "CSMD1", #*
    "FOXP1",
    "KIF21B",
    "PLK2",
    "PTK2B",
    "SCARA3",
    "SLC4A10", #*
    "SPG7"
)
fig5c_non_sig_genes <- c(
    "AP3B2",
    "ARHGAP44",
    "BBX",
    "BNIP3L",
    "CLU",
    "CNTN4",
    "CSDE1",
    "CXXC5",
    "DDHD2",
    "DGKZ",
    "MC1R",
    "FANCI",
    "HMOX2",
    "IGSF9B",
    "MAP3K11",
    "NDRG4",
    "RP11-677M14.7",
    "SOX2-OT",
    "SPATA33",
    "TBC1D5",
    "TXNRD1",
    "ZEB2"
)


## Get corresponding transcript IDs
sig_transcript_ids <- df_x_res_transcripts %>% 
    filter(gene_symbol %in% c(fig5c_sig_genes, fig5c_non_sig_genes)) %>% 
    pull(ensembl_transcript_id)


## Ensembl query
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl",
                      version = 107 # to align with GTEx v10
)
transcript_data <- getBM(
    attributes = c(
        "ensembl_transcript_id",
        "ensembl_transcript_id_version",
        "external_transcript_name",
        "ensembl_gene_id",
        "ensembl_gene_id_version",
        "external_gene_name",
        "transcript_biotype",
        "transcript_length",
        "description",
        "transcript_appris",
        "transcript_gencode_basic",
        #"transcript_primary_basic",
        "transcript_is_canonical",
        "transcript_mane_select",
        "transcript_mane_plus_clinical",
        "ccds"
    ),
    filters = "ensembl_transcript_id",
    values = sig_transcript_ids,
    mart = ensembl
) %>% 
    as_tibble() %>% 
    dplyr::rename("gene_symbol" = "external_gene_name", "transcript_symbol" = "external_transcript_name")
transcript_data


## Exon data
exon_data <- getBM(
    attributes = c(
        "ensembl_gene_id",
        "external_gene_name",
        "ensembl_transcript_id",
        #"external_transcript_name",
        "ensembl_exon_id",
        "chromosome_name",
        "strand",
        "exon_chrom_start",
        "exon_chrom_end",
        "rank"
    ),
    filters = "ensembl_transcript_id",
    values = sig_transcript_ids,
    mart = ensembl
) %>% 
    as_tibble() %>% 
    dplyr::rename("gene_symbol" = "external_gene_name",
                  "exon_start" = "exon_chrom_start",
                  "exon_end" = "exon_chrom_end"
    )
exon_data


## Regulatory regions
reg_mart <- useEnsembl(biomart = "ENSEMBL_MART_FUNCGEN", dataset = "hsapiens_regulatory_feature")
# Example: get regulatory regions for chromosome 8 from 3.35Mb to 4.1Mb
reg_regions <- getBM(
    attributes = c(
        "regulatory_stable_id", "feature_type_name", "feature_type_description",
        "chromosome_name", "chromosome_start", "chromosome_end", 
        "bound_seq_region_start", "bound_seq_region_end", 
        "epigenome_name", "epigenome_description"
    ),
    filters = c("chromosome_name", "start", "end"),
    values = list("8", 2832372, 5097953),
    mart = reg_mart
) %>% 
    as_tibble()


## GTEX
sig_gene_id_version <- transcript_data %>% 
    pull(ensembl_gene_id_version) %>% 
    unique
df_transcript_expression_gtex <- get_median_transcript_expression(
    gencodeIds = sig_gene_id_version,
    datasetId = "gtex_v10",
    tissueSiteDetailIds = "Brain_Anterior_cingulate_cortex_BA24",
    itemsPerPage = 100000
) %>% 
    as_tibble() %>% 
    clean_names() %>%
    mutate(ensembl_transcript_id = str_remove(transcript_id, "\\..*")) %>% 
    dplyr::rename("ensembl_transcript_id_version" = "transcript_id") %>% 
    
    left_join(transcript_data %>% dplyr::select(ensembl_transcript_id, ensembl_transcript_id_version, transcript_symbol, transcript_biotype)) %>%
    distinct()

## xQTL data
isoQTL_data <- read_table(paste0(base_dir, "data/DER-10a_hg19_isoQTL.significant.txt"))
df_isoQTL <- isoQTL_data %>% 
    dplyr::rename("ensembl_transcript_id_version" = "transcript_id", "position" = "SNP_start") %>%
    mutate(chromosome = str_remove(gene_chr, "chr") %>% as.numeric) %>% 
    dplyr::select(ensembl_transcript_id_version, chromosome, position, FDR) %>% 
    left_join(transcript_data %>% dplyr::select(ensembl_transcript_id_version, transcript_symbol, gene_symbol)) %>% 
    filter(!is.na(gene_symbol))

eQTL_data <- read_table(paste0(base_dir, "data/DER-08a_hg19_eQTL.significant.txt"))
df_eQTL <- eQTL_data %>% 
    dplyr::rename("ensembl_gene_id_version" = "gene_id", "position" = "SNP_start") %>%
    mutate(chromosome = str_remove(gene_chr, "chr") %>% as.numeric) %>%
    mutate(ensembl_gene_id = str_remove(ensembl_gene_id_version, "\\..*")) %>% 
    dplyr::select(ensembl_gene_id, chromosome, position, FDR) %>% 
    left_join(transcript_data %>% dplyr::select(ensembl_gene_id, transcript_symbol, gene_symbol)) %>% 
    filter(!is.na(gene_symbol))



# Write function to plot --------------------------------------------------


f_plot_isoform_data <- function(gene_of_interest, 
                                scaling_value = 1e6,
                                x_min = NULL, 
                                x_max = NULL
                                ) {
    
    # Subset and preprocess data for the current gene
    transcript_order <- df_x_res_transcripts %>% 
        filter(gene_symbol == gene_of_interest) %>% 
        arrange(pearsons_r) %>% 
        pull(transcript_symbol) %>% 
        unique()
    
    # Filter exons for gene of interest and calculate introns
    plot_data <- exon_data %>% 
        filter(gene_symbol == gene_of_interest) %>% 
        left_join(transcript_data) %>% 
        mutate(transcript_symbol = factor(transcript_symbol, levels = transcript_order)) %>% 
        group_by(transcript_symbol) %>%
        arrange(rank, .by_group = TRUE) %>%
        mutate(
            next_exon_start = lead(exon_start),
            intron_start = exon_end,
            intron_end = lead(exon_start),
            has_intron = !is.na(next_exon_start)
        ) %>%
        ungroup() %>% 
        mutate(is_principal = ifelse(str_detect(transcript_appris, "principal"), "*", ""))
    
    # Determine minimum and maximum values of plots
    chromosome_positions <- c(
        plot_data$intron_start,
        plot_data$intron_end,
        plot_data$exon_start,
        plot_data$exon_end,
        df_scz_broad_variants %>% filter(gene_symbol == gene_of_interest) %>% pull(position),
        df_eQTL %>% filter(gene_symbol == gene_of_interest) %>% pull(position)
    )
    chromosome_positions <- chromosome_positions[!is.na(chromosome_positions) & is.finite(chromosome_positions)]
    scale_min <- if (is.null(x_min)) min(chromosome_positions) / scaling_value else x_min
    scale_max <- if (is.null(x_max)) max(chromosome_positions) / scaling_value else x_max
    
    # Add arrow for direction of transcription
    arrow_data <- tibble(
        # start = ifelse(plot_data$strand > 0, min(plot_data$exon_start)/scaling_value, max(plot_data$exon_end)/scaling_value),
        # stop = ifelse(plot_data$strand > 0, max(plot_data$exon_end)/scaling_value, min(plot_data$exon_start)/scaling_value),
        start = ifelse(plot_data$strand > 0, scale_min, scale_max),
        stop = ifelse(plot_data$strand > 0, scale_max, scale_min),
        y = 0
    ) %>% distinct()
    
    # Exon plot
    p_exons <- plot_data %>% 
        mutate(transcript_biotype = str_replace_all(transcript_biotype, "_", " ")) %>% 
        ggplot(aes(y = transcript_symbol)) +
        geom_segment(aes(x = intron_start/scaling_value, xend = intron_end/scaling_value, y = transcript_symbol), linewidth = 0.25) +
        geom_rect(aes(
            xmin = exon_start/scaling_value, xmax = exon_end/scaling_value,
            ymin = as.numeric(transcript_symbol) - 0.4, ymax = as.numeric(transcript_symbol) + 0.4,
            fill = transcript_biotype), color = "black", linewidth = 0.1) +
        geom_arrowsegment(data = arrow_data, 
                          aes(x = start, xend = stop, y = y, yend = y),
                          arrow_positions = seq(0.01, 1.0, 0.09),
                          arrows = arrow(length = unit(0.25, "cm")),
                          color = "black", linewidth = 0.5) +
        geom_text(aes(x = -Inf, y = transcript_symbol, label = is_principal), hjust = -1, vjust = 0.75, size = 5) +
        labs(y = NULL, x = "Chromosome position (Mb)") +
        xlim(c(scale_min, scale_max)) +
        coord_cartesian(clip = "off") +
        theme(legend.position = "bottom", 
              legend.title = element_blank(), 
              legend.justification = "center",
              legend.key.size = unit(0.25, "cm"),
              legend.margin = margin(t = -10)
        )
    
    # Expression barplot
    p_expression <- df_transcript_expression_gtex %>% 
        filter(gene_symbol == gene_of_interest) %>% 
        mutate(transcript_symbol = factor(transcript_symbol, levels = transcript_order)) %>% 
        filter(!is.na(transcript_symbol)) %>% 
        ggplot(aes(x = median, y = transcript_symbol)) +
        geom_col(aes(fill = transcript_biotype), color = "black", linewidth = 0.25) +
        labs(x = "GTEx TPM", y = NULL) +
        coord_cartesian(clip = "off") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              legend.position = "none", axis.line = element_blank())
    
    # GWAS variants
    p_variants <- df_scz_broad_variants %>% 
        filter(gene_symbol == gene_of_interest) %>% 
        ggplot() +
        geom_point(aes(x = position/scaling_value, y = -log10(pval), fill = -log10(pval + 0.0001)),
                   color = "black", shape = 21, size = 1.5) +
        scale_fill_gradientn(colors = rev(brewer.rdylgn(100)), guide = "none") +
        labs(x = NULL, y = "-log10(p-value)", title = "GWAS loci") +
        xlim(c(scale_min, scale_max)) +
        coord_cartesian(clip = "off") +
        theme(plot.title = element_text(face = "plain"),
              axis.text.x = element_blank())
    
    # # Regulatory elements (brain only) ## take this out for now
    # brain_keywords <- c("brain", "astrocyte", "neural", "neuron", "nerve", "glia", "SK-N-MC")
    # p_reg <- reg_regions %>%
    #     filter(str_detect(epigenome_name, paste0(brain_keywords, collapse = "|"))) %>% 
    #     mutate(feature_type_name = factor(feature_type_name)) %>% 
    #     ggplot(aes(y = feature_type_name)) +
    #     geom_rect(aes(xmin = chromosome_start/scaling_value, xmax = chromosome_end/scaling_value, 
    #                   ymin = as.numeric(feature_type_name) - 0.4,
    #                   ymax = as.numeric(feature_type_name) + 0.4)) +
    #     labs(y = NULL, title = "Regulatory regions") +
    #     xlim(c(scale_min, scale_max)) +
    #     coord_cartesian(clip = "off") +
    #     theme(plot.title = element_text(face = "plain"))
    
    # eQTLs
    p_eQTL <- df_eQTL %>% 
        filter(gene_symbol == gene_of_interest) %>%
        dplyr::select(position, FDR) %>% distinct() %>% 
        
        ggplot() +
        geom_point(
            aes(x = position/scaling_value, y = -log10(FDR), fill = -log10(FDR + 0.0001)),
            color = "black", shape = 21, size = 1.5
        ) +
        scale_fill_gradientn(colors = rev(brewer.rdylgn(100)), guide = "none") +
        labs(x = NULL, title = "eQTLs") +
        xlim(c(scale_min, scale_max)) +
        coord_cartesian(clip = "off") +
        theme(plot.title = element_text(face = "plain"),
              axis.text.x = element_blank())
    
    
    ## Combine & plot
    layout <- c(
        patchwork::area(t = 1, b = 90, l = 1, r = 140),
        patchwork::area(t = 91, b = 180, l = 1, r = 140),
        patchwork::area(t = 181, b = 450, l = 1, r = 140),
        patchwork::area(t = 181, b = 450, l = 141, r = 180)
    )
    
    # Combine and print plot
    full_plot <- p_variants + p_eQTL + p_exons + p_expression +
        patchwork::plot_layout(design = layout) +
        patchwork::plot_annotation(title = paste0(gene_of_interest, " (chr ", unique(plot_data$chromosome_name), ")"))
    
    print(full_plot)
    
}


f_plot_isoform_data("CSMD1")
f_plot_isoform_data("ARHGAP44")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "D.ARHGAP44", .x), width = 6.5, height = 4.5)
)

f_plot_isoform_data("PTK2B")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "D.PTK2B", .x), width = 6.5, height = 4.5)
)

f_plot_isoform_data("CSDE1")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "D.CSDE1", .x), width = 6.5, height = 4.5)
)




# Loop through and write PDF of all plots ---------------------------------


# Define your vector of gene symbols
genes_of_interest <- c(fig5c_sig_genes, fig5c_non_sig_genes)

# Scale factor (to divide x-axis)
scaling_value <- 1e6

# Open a single PDF device
pdf("~/Desktop/gene_isoform_plots.pdf", width = 12, height = 10)

for (gene_of_interest in genes_of_interest) {
    
    print(gene_of_interest)
    f_plot_isoform_data(genes_of_interest) %>% print
    
}

# Close the PDF device
dev.off()



# Save final plot for markdown --------------------------------------------



fig5c <- f_plot_isoform_data("CSMD1")
save(fig5c, file = paste0(analysis_objects_dir, "Fig5c.Rdata"))

figS14a <- f_plot_isoform_data("PTK2B")
figS14b <- f_plot_isoform_data("ARHGAP44")
save(figS14a, figS14b, file = paste0(analysis_objects_dir, "FigS14.RDS"))


