
################################################################################

# Put together figures for Fig 1 analysis pipeline

################################################################################


# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "08Mar2024_GENES_qSVAgeSexRaceGC"

## DEFINE WGCNA PARAMETERS
soft_power <- 3
minimum_size <- 40
tree_cut_height <- 0.98

# LOAD OBJECTS 
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules
load(paste0(base_dir, "objects/hsapiens_genome_v110.RDS")) # df_hsapiens_genome

# IDENTIFY MODULE SET OF INTEREST  
df_modules_filt <- df_modules %>% 
    filter(min_size == minimum_size & cut_height == tree_cut_height) %>% 
    unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module)


# WGCNA heatmap & dendrogram -------------------------------------------------------------------


# CONVERT REGRESSED VSD COUNTS TO TRANSPOSED MATRIX
m_vsd_regress <- df_vsd_regress %>% column_to_rownames("sample")

# CALCULATE CO-EXPRESSION SIMILARITY AND ADJACENCY
set.seed(20240306)
adjacency <- adjacency(m_vsd_regress, type = "signed hybrid", power = soft_power)

# TOPOLOGICAL OVERLAP MATRIX (TOM)
set.seed(20240306)
TOM <- TOMsimilarity(adjacency, TOMType = "signed")

# subset for plotting
baby_TOM <- TOM[1:2000, 1:2000]
baby_diss_TOM <- 1 - baby_TOM

# HIERARCHICAL CLUSTERING ON TOM  
baby_gene_tree <- hclust(as.dist(baby_diss_TOM))

# plot the resulting clustering tree
png(file = paste0(base_dir, "outputs/figures/for_manuscript/1.WGCNA_dendrogram.png"), width = 1200, height = 700)
par(mfrow = c(1, 1))
plot(baby_gene_tree, labels = FALSE, axes = FALSE, xlab = "", ylab = "", sub = "", main = "")
dev.off()

# PLOT NETWORK HEATMAP
colnames(TOM) <- colnames(adjacency)
rownames(TOM) <- rownames(adjacency)
diag(TOM) <- NA
annotation <- df_modules_filt %>% 
    dplyr::select(ensembl_gene_id, module) %>% 
    filter(ensembl_gene_id %in% colnames(baby_TOM)) %>% 
    mutate(module = as.character(module)) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = colnames(baby_TOM))) %>% 
    arrange(ensembl_gene_id) %>% 
    column_to_rownames("ensembl_gene_id")
all(rownames(annotation) == rownames(baby_TOM))
module_colors <- df_modules_filt %>% 
    mutate(module = as.character(module)) %>% 
    filter(ensembl_gene_id %in% colnames(baby_TOM)) %>% 
    dplyr::select(module, color) %>% 
    distinct %>% 
    deframe
annotation_col = list(module = module_colors)
pheatmap(
    mat = baby_TOM,
    color = colorRampPalette(c("white", "black"))(100),
    #color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_row = 0,
    treeheight_col = 0,
    annotation_row = annotation,
    annotation_col = annotation,
    annotation_colors = annotation_col,
    annotation_names_row = FALSE, 
    annotation_names_col = FALSE,
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = FALSE,
    legend = FALSE,
    annotation_legend = FALSE,
    filename = paste0(base_dir, "outputs/figures/for_manuscript/1.WGCNA_heatmap.png"),
    width = 2,
    height = 1.75
)




# Differential expression -------------------------------------------------


tibble(
    dx = "Control",
    expression = rnorm(n = 50, mean = -1, sd = 1)
) %>% 
    bind_rows(
        tibble(
            dx = "SCZ",
            expression = rnorm(n = 50, mean = 1, sd = 1)
        )
    ) %>% 
    
    ggplot(aes(x = expression, y = dx, color = dx)) +
    geom_point(position = position_jitter(width = 0.05), shape = 21, size = 1.5) +
    geom_boxplot(fill = "transparent", outlier.shape = NA) +
    scale_color_manual(values = c("Control" = "#9c522c", "SCZ" = "#d07444"), guide = "none") +
    labs(x = "", y = "", title = "Differential expression") +
    theme_classic() +
    theme(axis.text.x = element_blank())

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/1.DE", .x),
                  width = 3, height = 2)
)



# Alternative splicing ----------------------------------------------------


## IDENTIFY BIOTYPE FOR EACH TRANSCRIPT ID
df_transcript_biotype <- df_hsapiens_genome %>% 
    filter(type == "transcript") %>% 
    dplyr::select(transcript_id, transcript_biotype)

# IDENTIFY EXONS FOR SIGNIFICANT TRANSCRIPTS
exons <- df_hsapiens_genome %>% 
    filter(gene_name == "CHGA") %>% 
    mutate(transcript_name = ifelse(is.na(transcript_name), transcript_id, transcript_name)) %>% 
    filter(!is.na(transcript_name)) %>% 
    filter(type == "exon") %>% 
    mutate(start = start/1e6,
           end = end/1e6) %>% 
    mutate(transcript_name = str_remove(transcript_name, ".*-"),
           transcript_biotype = str_replace_all(transcript_biotype, "_", " ")
    )

# example plot (main text)
exons %>% 
    
    # plot isoforms
    ggplot(aes(
        xstart = start,
        xend = end,
        y = transcript_name
    )) +
    geom_range(fill = "#4C768A", color = "black") +
    geom_intron(
        data = to_intron(exons %>% filter(gene_name == "CHGA"), "transcript_name"),
        aes(strand = strand)
    ) +
    #facet_wrap(vars(gene_name)) +
    guides(fill = guide_legend(title = "Biotype", title.position = "top", title.hjust = 0.5)) +
    labs(y = "Isoform ", x = "Position (Mb)"
    ) +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank()
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/1.AS", .x),
                  width = 5, height = 2)
)



# sample network ----------------------------------------------------------


# random graph
net <- rgraph(50, mode = "graph", tprob = 0.2)
net <- network(net, directed = FALSE)
ggnet2(net,
       node.color = "#d07444",
       node.size = 2.5,
       edge.color = "#d07444",
       edge.alpha = 0.5
) +
    theme_void()

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/1.gene_net", .x),
                  width = 3, height = 2)
)

ggnet2(net,
       node.color = "#4C768A",
       node.size = 2.5,
       edge.color = "#4C768A",
       edge.alpha = 0.5
) +
    theme_void()

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/1.transcript_net", .x),
                  width = 3, height = 2)
)



# gene weight vector ------------------------------------------------------


## LOAD DATA

# genes
df_x_res_genes <- read_excel(paste0(base_dir, "outputs/tables/for_manuscript/TableS3_GENES_GRCCA_res_withGray.xlsx"), sheet = 3) %>% 
    mutate(module = factor(module, levels = 0:max(as.numeric(module))))

# transcripts
df_y_res_transcripts <- read_excel(paste0(base_dir, "outputs/tables/for_manuscript/TableS7_TRANSCRIPTS_GRCCA_res_withGray.xlsx"), sheet = 2)
df_x_res_transcripts <- read_xlsx(paste0(base_dir, "outputs/tables/for_manuscript/TableS7_TRANSCRIPTS_GRCCA_res_withGray.xlsx"), sheet = 3) %>% 
    mutate(module = factor(module, levels = 0:max(as.numeric(module))))

## PLOT GENES
df_x_res_genes %>% 
    
    ggplot(aes(x = pearsons_r, y = reorder(ensembl_gene_id, pearsons_r))) +
    geom_col(aes(fill = pearsons_r)) +
    geom_vline(xintercept = 0, lty = 2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
    guides(fill = guide_colorbar(title = "Pearson's r")) +
    #scale_color_manual(values = c("*" = "black", "-" = "transparent")) +
    coord_flip() +
    labs(y = "Gene", x = "Pearson's r") +
    theme_classic() +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(),
          legend.position = "none"
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/1.gene_vector", .x),
                  width = 5, height = 1.1)
)

## PLOT GENES
df_x_res_transcripts %>% 
    
    ggplot(aes(x = pearsons_r, y = reorder(ensembl_transcript_id, pearsons_r))) +
    geom_col(aes(fill = pearsons_r)) +
    geom_vline(xintercept = 0, lty = 2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
    guides(fill = guide_colorbar(title = "Pearson's r")) +
    #scale_color_manual(values = c("*" = "black", "-" = "transparent")) +
    coord_flip() +
    labs(y = "Isoform", x = "Pearson's r") +
    theme_classic() +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(),
          legend.position = "none"
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/1.transcript_vector", .x),
                  width = 5, height = 1.1)
)



# module enrichment -------------------------------------------------------

df_modules_count <- df_modules_filt %>% 
    count(module, color)
module_colors <- df_modules_count %>% 
    dplyr::select(module, color) %>% 
    deframe

# random graph
l_graphs <- list()
for (mod in 1:10) {
    
    mod_color <- module_colors[paste0(mod)]
    n_nodes <- floor((df_modules_count %>% filter(module == mod) %>% pull(n))/50)
    net <- rgraph(n_nodes, mode = "graph", tprob = 0.25)
    net <- network(net, directed = FALSE)
    p_graph <- ggnet2(net,
           node.color = mod_color,
           node.size = 2.5,
           edge.color = mod_color,
           edge.alpha = 0.6
    ) +
        theme_void()
    
    l_graphs[[mod]] <- p_graph
    
    # print & save
    p_graph
    map(
        .x = c(".png", ".pdf"),
        .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/1.gene_mod_nets_M", mod, .x),
                      width = 1.5, height = 1)
    )
    
}
wrap_plots(l_graphs)






# gene/transcript enrichment ---------------------------------------------------------

noise <- 2000
tibble(
    x = jitter(seq(1:max_n), noise),
    y = jitter(x^2, noise)
) %>% 
    bind_rows(
        tibble(
            x = jitter(seq(1:max_n), noise),
            y = jitter(-x^2, noise)
        )
    ) %>% 
    filter(y > 0) %>% 
    ggplot(aes(x = y, y = x)) +
    geom_point(fill = "#d07444", shape = 21, aes(alpha = y, size = y)) +
    labs(x = "-log10(p)", y = "NES") +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/1.gene_enrich", .x),
                  width = 3, height = 2)
)

## TRANSCRIPTS
tibble(
    x = jitter(seq(1:max_n), noise),
    y = jitter(x^2, noise)
) %>% 
    bind_rows(
        tibble(
            x = jitter(seq(1:max_n), noise),
            y = jitter(-x^2, noise)
        )
    ) %>% 
    filter(y > 0) %>% 
    ggplot(aes(x = y, y = x)) +
    geom_point(fill = "#4C768A", shape = 21, aes(alpha = y, size = y)) +
    labs(x = "-log10(p)", y = "NES") +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/1.transcript_enrich", .x),
                  width = 3, height = 2)
)



