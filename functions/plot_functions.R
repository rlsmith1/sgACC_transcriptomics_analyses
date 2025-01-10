
###################################################################

# Functions to generate all plots from these analyses

###################################################################


# Plot theme --------------------------------------------------------------

## Color themes
dx_colors <- c("Control" = "#0072B2", "BD" = "#9966FF", "MDD" = "#009E73", "SCZ" = "#E69F00")
ontology_colors <- c("BP" = "#F8766D", "CC" = "#00BA38", "MF" = "#619CFF")
gene_weight_color_scale <- brewer.rdbu(100) %>% rev

## Plot aesthetics
theme_set(
    theme_cowplot() +
        theme(plot.title = element_text(size = 11),
              strip.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 9),
              legend.title = element_text(size = 9),
              legend.text = element_text(size = 8),
              legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margin
              plot.margin = margin(t = 0, r = 5, b = 0, l = 0)
        )
)


# plot GO results by ontology ---------------------------------------------

f_plot_go_by_ontology <- function(df_go_res,
                                  n_top_pathways = 10,
                                  pathway_text_width = 40,
                                  pathway_text_size = 9,
                                  plot_title = "GO pathway results",
                                  n_facet_rows = 2,
                                  legend_position = c(0.7, 0.25)
) {
    df_go_res %>% 
        
        # clean dataframe (if not already)
        clean_names %>% 
        
        # take top 10 paths by p-value per ontology to plot
        group_by(ontology) %>% 
        arrange(-p_adjust) %>% 
        mutate(row = row_number(),
               gene_ratio = parse(text = gene_ratio) %>% eval
        ) %>% 
        top_n(n = n_top_pathways, wt = row) %>% 
        
        # indicate whether pathway survived FDR correction
        mutate(`survives FDR` = ifelse(p_adjust < 0.05, "yes", "no")) %>% 
        mutate(description = str_wrap(description, width = pathway_text_width)) %>%
        
        # plot
        ggplot(aes(x = -log10(p_adjust), y = reorder_within(description, within = ontology, by = -p_adjust))) +
        geom_point(aes(size = gene_ratio, fill = ontology, color = `survives FDR`),
                   shape = 21, stroke = 1) +
        geom_vline(aes(xintercept = -log10(0.05)), color = "gray") +
        facet_wrap(vars(ontology), scales = "free", nrow = n_facet_rows) +
        scale_size_continuous(range = c(2, 4)) +
        scale_color_manual(values = c("yes" = "black", "no" = "transparent"), guide = "none") +
        scale_fill_manual(values = ontology_colors, guide = "none") +
        scale_y_reordered() +
        coord_cartesian(clip = "off") +
        labs(y = "", x = "log10(FDR)",
             title = plot_title
        ) +
        guides(size = guide_legend(title = "Gene ratio")) +
        theme(
            axis.text.y = element_text(size = pathway_text_size),
            legend.position = legend_position
        )
}


# plot GSEA results by ontology -------------------------------------------

f_plot_gsea_by_ontology <- function(df_gsea_res,
                                    ontologies = c("BP", "MF", "CC"),
                                    n_pathways_labeled = 5,
                                    pathway_text_size = 3,
                                    pathway_text_width = 30,
                                    pathway_text_overlaps = 40,
                                    pathway_box_padding = 0.75,
                                    orientation = c("vertical", "horizontal"),
                                    plot_title = "GSEA enrichments",
                                    legend_position = c(0.87, 0.97)
) {
    if (orientation == "vertical") {ncol <- 1} else if (orientation == "horizontal") {ncol <- 3} else (stop("Please select orientation for panels"))
    
    # Generate GSEA plot
    p <- df_gsea_res %>% 
        dplyr::filter(ontology %in% ontologies) %>%
        mutate(enrichment_sign = ifelse(nes < 0, "neg", "pos")) %>% 
        mutate(nes_fill = ifelse(p_adjust > 0.05, NA_real_, nes)) %>% 
        group_by(ontology, enrichment_sign) %>% 
        mutate(order = row_number(),
               label = ifelse(order %in% seq(1, n_pathways_labeled) | -log10(p_adjust) > 10, description %>% str_wrap(pathway_text_width), NA)
        ) %>% 
        
        # plot
        ggplot(aes(x = -log10(p_adjust), y = nes)) +
        geom_point(aes(size = set_size, color = nes_fill, alpha = abs(nes))) +
        geom_vline(aes(xintercept = -log10(0.05)), color = "black", lty = 2) +
        geom_hline(aes(yintercept = 0), color = "gray") +
        geom_text_repel(aes(label = label), min.segment.length = 0, size = pathway_text_size,
                        box.padding = pathway_box_padding, max.overlaps = pathway_text_overlaps) +
        facet_wrap(vars(ontology), ncol = ncol) +
        scale_size_continuous(range = c(0.5, 4), limits = c(10, 500)) +
        scale_alpha_continuous(range = c(0.1, 1), guide = "none") +
        scale_color_gradientn(colors = rev(brewer.rdbu(100)), na.value = "gray", guide = "none") +
        guides(size = guide_legend(title = "n genes")) +
        labs(y = "Normalized enrichment score", x = "-log10(FDR)",
             title = plot_title
        ) +
        theme(legend.position = legend_position,
              legend.key.size = unit(0.2, "cm")
        )
    
    # Remove strip text if there's only one ontology
    if (length(ontologies) == 1) {p <- p + theme(strip.text = element_blank())}
    
    # Return plot
    return(p)
}



# Plot gene subset position in z-score distribution -----------------------

f_plot_gene_position <- function(df, 
                                 df_subset,
                                 x_axis_variable = c("z-score", "structure correlation", "l2fc"),
                                 plot_title = "Gene position in distribution",
                                 gene_size = 2.5, 
                                 max_text_overlaps = 20,
                                 gene_box_padding = 0.75
) {
    
    # Determine which variable to plot on the x-axis
    if (x_axis_variable == "z-score") {
        df <- df %>% dplyr::rename("x" = "z_score")
        df_subset <- df_subset %>% dplyr::rename("x" = "z_score")
        xlab <- "Gene weight z-score"
        vline_intercepts <- c(-2, 2)
    } else if (x_axis_variable == "structure correlation") {
        df <- df %>% dplyr::rename("x" = "pearsons_r")
        df_subset <- df_subset %>% dplyr::rename("x" = "pearsons_r")
        xlab <- "Gene structure correlation"
        cor_thresh <- df[which.min(abs(0.05 - df$p_adj)), ] %>% pull(x) %>% abs
        vline_intercepts <- c(-cor_thresh, cor_thresh)
    } else if (x_axis_variable == "l2fc") {
        df <- df %>% dplyr::rename("x" = "log2fold_change")
        df_subset <- df_subset %>% dplyr::rename("x" = "log2fold_change")
        xlab <- "Gene L2FC"
        vline_intercepts <- 0
    }
    
    # Set scale max and min depending on variable being plotted
    scale_max <- max(df %>% pull(x) %>% min, df %>% pull(x) %>% max) %>% abs
    
    # Plot
    df %>% 
        ggplot(aes(x = x)) +
        geom_vline(xintercept = vline_intercepts, color = "gray") +
        geom_density() +
        geom_point(data = df_subset %>% arrange(significant), shape = 21, size = gene_size,
                   aes(x = x, y = 0, fill = x, color = as.factor(significant))) +
        geom_text_repel(data = df_subset %>% dplyr::filter(significant == 1), 
                        min.segment.length = 0.1, box.padding = gene_box_padding, size = 2.5, max.overlaps = max_text_overlaps,
                        aes(x = x, y = 0, label = gene_symbol)
        ) +
        scale_fill_gradientn(colors = gene_weight_color_scale, 
                             values = rescale(c(-scale_max, 0, scale_max)),
                             limits = c(-scale_max, scale_max),
                             na.value = "lightgray", guide = "none") +
        scale_color_manual(values = c("transparent", "yellow"), guide = "none") +
        coord_flip() +
        xlim(c(-scale_max, scale_max)) +
        labs(x = xlab, y = "", title = plot_title)
    
}



# Add genes to an existing density distribution plot ----------------------

f_add_genes <- function(df, 
                        list_no,
                        #jitter_position = position_jitter(width = 0.1, height = 0, seed = 123),
                        n_genes_to_label = 10,
                        label_x_pos = 0.1,
                        point_size = 3,
                        y_axis_scaling = -10
                        ) {
    
    return(
        list(
            geom_point(
                data = df %>% 
                    dplyr::filter(ensembl_gene_id %in% benchmarking_lists[[list_no]]) %>% 
                    mutate(significant = pvalue < 0.05, 1, 0) %>% 
                    arrange(significant, x),
                mapping = aes(x = x, y = list_no*y_axis_scaling, fill = x, color = significant), 
                position = position_jitter(width = 0.1, seed = 123), 
                shape = 21, size = point_size
            ),
            geom_text_repel(data = df %>% 
                                dplyr::filter(ensembl_gene_id %in% benchmarking_lists[[list_no]]) %>%  #& pvalue < 0.05) %>%
                                top_n(n = n_genes_to_label/2, wt = abs(x)),
                            mapping = aes(x = x, y = list_no*y_axis_scaling, label = gene_symbol),
                            position = position_jitter(width = 0.1, seed = 123), 
                            size = 3, box.padding = 0.75, max.overlaps = 20, force = 10, min.segment.length = 0
            ),
            annotate(geom = "text", x = label_x_pos, y = list_no*y_axis_scaling, hjust = 0, fontface = 2, size = 3.25,
                     label = names(benchmarking_lists[list_no]) %>% str_replace_all("_", " ")
            )
        )
    )
    
}



# Plot module hypergeometric test results ---------------------------------

f_plot_module_hypergeometric <- function(df_hypergeometric_res,
                                         x_axis = c("p-value", "p-adj"),
                                         plot_title = "Module overrepresentation",
                                         #subset_x_labs = TRUE,
                                         include_module_labels = TRUE,
                                         module_text_size = 3,
                                         legend_position = "none",
                                         p_value_threshold = 0.05,
                                         include_threshold_label = TRUE
) {
    
    # specify whether to plot the unadjusted or adjusted p-value on the y-axis
    # if there are no genes in the overlap, set the p-value to 1
    if (x_axis == "p-value") {
        df_hypergeometric_res <- df_hypergeometric_res %>% dplyr::rename("p" = "p_value") %>% mutate(p = ifelse(overlap_n == 0, 1, p))
        x_axis_title = "-log10(Hypergeometric p-value)"
        threshold_label = paste0("p < ", p_value_threshold)
    } else if (x_axis == "p-adj") {
        df_hypergeometric_res <- df_hypergeometric_res %>% dplyr::rename("p" = "p_adj") %>% mutate(p = ifelse(overlap_n == 0, 1, p))
        x_axis_title = "-log10(Hypergeometric FDR)"
        threshold_label = paste0("FDR < ", p_value_threshold)
    } else(stop("Must plot either p-value or p-adj"))
    
    # include text at the horizontal threshold line
    if (include_threshold_label == TRUE) {
        threshold_label <- threshold_label
    } else if (include_threshold_label == FALSE) {
        threshold_label <- ""
    }
    
    # plot
    p <- df_hypergeometric_res %>%
        ggplot(aes(x = -log10(p), y = module)) +
        geom_vline(aes(xintercept = -log10(p_value_threshold)), color = "gray", linewidth = 0.5) +
        geom_segment(aes(x = 0, xend = -log10(p), yend = module), linewidth = 0.25, color = "black") +
        geom_point(aes(fill = module, size = overlap_n/total_n, color = p < 0.05), shape = 21) +
        annotate(geom = "text", label = threshold_label, size = 3, hjust = -0.25,
                 y = 1, x = -log10(p_value_threshold)) +
        scale_fill_manual(values = module_colors, guide = "none") +
        scale_color_manual(values = c("transparent", "black")) +
        scale_y_discrete(limits = rev) +
        labs(x = x_axis_title, y = NULL,
             title = plot_title
        ) +
        coord_cartesian(clip = "off") +
        theme(legend.position = legend_position,
              legend.box = "horizontal",
        )
    
    # add module labels if specified
    if (include_module_labels == TRUE) {
        p <- p + geom_text(aes(label = module_label), hjust = 0.5, vjust = -1, size = module_text_size)
    }
    
    # Return plot
    return(p)
    
    # # subset x labels if specific
    # if (subset_x_labs == TRUE) {
    #     
    #     sig_modules <- df_hypergeometric_res %>% dplyr::filter(p < 0.05) %>% pull(module)
    #     x_labs <- map(
    #         .x = names(module_colors),
    #         .f = ~ if (.x %in% sig_modules) {.x} else {""}
    #     ) %>% unlist()
    #     p <- p + scale_x_discrete(labels = x_labs)
    #     return(p)
    #     
    # } else {return(p)}
    
}

