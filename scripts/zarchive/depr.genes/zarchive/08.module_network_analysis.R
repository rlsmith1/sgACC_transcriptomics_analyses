

########################################################################################

# Characterize network properties of control modules in diagnostic groups

########################################################################################


# libraries ---------------------------------------------------------------


  library(tidyverse)
  library(WGCNA)
  library(igraph)
  library(patchwork)
  library(janitor)
  library(RColorBrewer)
  library(tidymodels)
  library(corrr)


# set theme for plots -----------------------------------------------------

  theme_set(theme_bw() +
              theme(plot.title = element_text(size = 18),
                    axis.title = element_text(size = 15),
                    axis.text = element_text(size = 12),
                    strip.text = element_text(size = 15),
                    legend.title = element_text(size = 15),
                    legend.text = element_text(size = 12)))
  



# data --------------------------------------------------------------------


  prefix <- "20221209_185samples_19kgenes_vsd_qSVA17_resids"
  ctrl_prefix <- "20221209_55samplesCTRL_19kgenes_vsd_qSVA17_resids"

# LOAD REGRESSED DATA
  load(paste0("objects/", prefix, ".RDS")) # df_vsd_regress

# FILTER FOR CONTROL ONLY
  df_control <- df_vsd_regress %>% 
    filter(grepl("control", sample))

# LOAD CONTROL MODULE ASSIGNMENTS
  load(paste0("objects/", ctrl_prefix, "_SFT4_MODS.RDS")) # df_modules

# IDENTIFY MODULE SET OF INTEREST  
  df_modules_filt <- df_modules %>% 
    filter(min_size == 50 & cut_height == 0.99) %>% 
    dplyr::select(ensembl_gene_id, module) %>% 
    mutate(module = factor(module, levels = 0:max(module)))

# COVARIATE DATA
  load("data/covariates/185_all_covariates_clean.Rdata")
  


# create graph based on gene correlations (not TOM weights) --------------

  
  df_mods_vsd <- df_modules_filt %>% 
    left_join(df_vsd_regress, by = "ensembl_gene_id") %>% 
    mutate(dx = gsub(".*_", "", sample) %>% 
             factor(levels = c("control", "bipolar", "mdd", "schizo")), 
           .before = 5)
  
# CALCULATE CORRELATIONS  IN EACH MODULE IN EACH DX
  df_mods_cors <- df_mods_vsd %>% 
    group_by(module, dx) %>% 
    nest() %>% 
    arrange(dx, module) %>% 
    filter(module != 0) %>% 
    
    mutate(cor = map(
      
      .x = data,
      .f = ~ .x %>% 
        
        pivot_wider(id_cols = sample,
                    names_from = ensembl_gene_id,
                    values_from = vsd) %>% 
        as.data.frame %>% 
        column_to_rownames("sample") %>% 
        cor %>% 
        as.data.frame() %>% 
        rownames_to_column("gene1") %>% 
        as_tibble %>% 
        pivot_longer(contains("ENSG"),
                     names_to = "gene2",
                     values_to = "r")
      
    )
    )
  
# CALCULATE CORRELATIONS IN EACH DX 
  df_dx_gene_cors <- df_vsd_regress %>%
    pivot_wider(id_cols = sample,
                names_from = ensembl_gene_id,
                values_from = vsd) %>% 
    
    mutate(dx = gsub(".*_", "", sample)  %>% 
             factor(levels = c("control", "bipolar", "mdd", "schizo")), 
           .before = 2) %>% 
    group_by(dx) %>% 
    nest() %>% 
    arrange(dx) %>% 
    
    mutate(cor = map(
      
      .x = data,
      .f = ~ .x %>% 
        
        as.data.frame %>% 
        column_to_rownames("sample") %>% 
        cor %>% 
        as.data.frame() %>% 
        rownames_to_column("gene1") %>% 
        as_tibble %>% 
        pivot_longer(contains("ENSG"),
                     names_to = "gene2",
                     values_to = "r")
      
    )
    )
  
  
# PLOT CORRELATIONS ACROSS DISORDERS  
  df_modules_filt
  df_dx_gene_cors %>% 
  
    ungroup %>% 
    filter(dx == "control") %>% 
    unnest(cols = c(cor)) %>% 
    filter(module == 14) %>% 
    
    group_by(module, dx, gene1) %>% 
    mutate(gene1_mean_r = mean(r)) %>% 
    arrange(module, dx, -gene1_mean_r, -r) %>% 
      
    mutate(gene1 = factor(gene1, levels = unique(.$gene1)),
           gene2 = factor(gene2, levels = unique(.$gene1))) %>% 
    mutate(r = ifelse(gene1 == gene2, NA_real_, r)) %>% 
    
    ggplot(aes(x = gene1, y = gene2, fill = r)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black", na.value = "black", limits = c(-1, 1)) +
    theme(axis.text.x = element_blank())
  
  
# CREATE GRAPH  
  df_mods_graphs <- df_mods_cors %>% 
    
    mutate(graph = map(
      
      .x = cor,
      .f = ~ .x %>% 
        filter(gene1 != gene2) %>% 
        graph_from_data_frame(directed = FALSE)
    )
    )
  

# calculate degree & hub scores -------------------------------------------

  
  df_mods_degree_hub <- df_mods_graphs %>% 
    mutate(degree = map(.x = graph,
                        .f = ~ graph.strength(.x)
    ),
    hub_score = map(.x = graph,
                    .f = ~ hub_score(.x)$vector
    )
    )  %>% 
    dplyr::select(-graph, -data) %>% 
    unnest(cols = c(degree, hub_score)) %>% 
    distinct()
    
  df_mods_degree_hub
  
  df_dx_graph %>% 
    ggplot(aes(x = hub_score)) +
    geom_density(aes(fill = dx), alpha = 0.5) +
    facet_wrap( ~ module, scales = "free") +
    ggtitle("Hub score distributions by module")
  
  
# calculate TOM for each diagnostic group to get gene weights -------------

  
# CALCULATE ADJACENCY AND TOM FOR EACH DIAGNOSTIC GROUP SEPARATELY (SFT = 4)  
  df_dx_wgcna <- df_vsd_regress %>% 
    pivot_wider(id_cols = sample, 
                names_from = ensembl_gene_id, 
                values_from = resids) %>% 
    mutate(dx = gsub(".*_", "", sample), .before = 1) %>% 
    group_by(dx) %>% 
    nest() %>% 
    
    mutate(data_mat = map(.x = data,
                          .f = ~ .x %>%     
                            as.data.frame %>% 
                            column_to_rownames("sample") %>% 
                            as.matrix)) %>% 

    mutate(adj = map(.x = data_mat,
                     .f = ~ adjacency(.x, power = 4)
                     )) %>% 
    
    mutate(tom = map(.x = adj,
                     .f = ~ TOMsimilarity(.x))
           ) %>% 
    dplyr::select(dx, tom)
  
  save(df_dx_wgcna,
       file = paste0("objects/", ctrl_prefix, "_4DX_SFT4_TOM.RDS"))
  
  

# PIVOT TOM INTO IGRAPH FORMAT
  load(paste0("objects/", ctrl_prefix, "_4DX_SFT4_TOM.RDS"))
  df_dx_tom <- df_dx_wgcna %>% 
    
    mutate(tom = map(.x = tom,
                     .f = ~ .x %>% 
                       as.data.frame %>% 
                       as_tibble %>% 
                       dplyr::rename_all( ~ unique(df_vsd_regress$ensembl_gene_id)) %>% 
                       mutate(gene1 = unique(df_vsd_regress$ensembl_gene_id), .before = 1) %>% 
                       pivot_longer(starts_with("ENSG"), names_to = "gene2", values_to = "weight") %>% 
                       
                       left_join(df_modules_filt %>% 
                                   dplyr::rename("gene1" = "ensembl_gene_id",
                                                 "module1" = "module"),
                                 by = "gene1"
                       ) %>% 
                       
                       left_join(df_modules_filt %>% 
                                   dplyr::rename("gene2" = "ensembl_gene_id",
                                                 "module2" = "module"),
                                 by = "gene2"
                       ) %>% 
                       
                       filter(module1 == module2) %>% 
                       dplyr::select(-module2) %>% 
                       dplyr::rename("module" = "module1")
    )
    
    )
  
  save(df_dx_tom,
       file = paste0("objects/", ctrl_prefix, "_4DX_SFT4_TOM.RDS"))
  
# PLOT TOM  
  df_dx_tom %>% 
    unnest(cols = c(tom)) %>% 
    
    # filter(dx == "bipolar") %>% pull(tom) %>% .[[1]] %>% 
    filter(gene1 != gene2) %>% 
    arrange(desc(weight)) %>% 
    mutate(thresh = ifelse(row_number() <= 0.1*nrow(.), 1, 0)) %>% 

    ggplot(aes(x = gene1, y = gene2, fill = thresh)) +
    geom_tile() +
    facet_wrap( ~ dx) + 
    scale_fill_gradient(low = "white", high = "black") +
    theme(axis.text = element_blank())
  
  
  tibble(x = seq(1, 10, 1)) %>% 
    arrange(desc(x)) %>% 
    mutate(y = ifelse(row_number() <= 0.1*nrow(.), 1, 0))

  
# create graph for each diagnostic group ----------------------------------

  
# CALCULATE MEDIAN CONNECTIVITY FOR EACH DIAGNOSIS
  df_dx_sft4 <- df_vsd_regress %>% 
    pivot_wider(id_cols = sample, 
                names_from = ensembl_gene_id, 
                values_from = resids) %>% 
    mutate(dx = gsub(".*_", "", sample), .before = 1) %>% 
    group_by(dx) %>% 
    nest() %>% 
    
    mutate(data_mat = map(.x = data,
                          .f = ~ .x %>%     
                            as.data.frame %>% 
                            column_to_rownames("sample") %>% 
                            as.matrix)) %>% 
    mutate(sft = map(.x = data,
                     .f = ~ pickSoftThreshold(.x, powerVector = 4, verbose = 5))) %>% 
    mutate(median_k = map(.x = sft,
                          .f = ~ .x$fitIndices$median.k),
           r_sq = map(.x = sft,
                          .f = ~ .x$fitIndices$SFT.R.sq)) %>% 
    unnest(cols = c(median_k,r_sq))
    

# GENERATE GRAPH FOR EACH MODULE IN EACH DX GROUP  
  load(paste0("objects/", ctrl_prefix, "_4DX_SFT4_TOM.RDS"))
  
  df_dx_graph <- df_dx_tom %>% 
    unnest(cols = c(tom)) %>% 
    group_by(dx, module) %>% 
    nest() %>% 
    mutate(graph = map(.x = data,
                       .f = ~ graph_from_data_frame(.x, directed = FALSE))) %>% 
    dplyr::select(dx, module, graph) %>% 
    mutate(dx = factor(dx, levels = c("control", "bipolar", "mdd", "schizo"))) %>% 
    arrange(dx, module)
  

# CALCULATE DEGREE AND HUB SCORE
  df_dx_graph <- df_dx_graph %>% 
    mutate(degree = map(.x = graph,
                        .f = ~ graph.strength(.x)
    ),
    hub_score = map(.x = graph,
                    .f = ~ hub_score(.x)$vector
    )
    )  %>% 
    unnest(cols = c(degree, hub_score))
  
  save(df_dx_graph,
       file = paste0("objects/", ctrl_prefix, "_4DX_SFT4_GRAPH.RDS"))
  
# NORMALIZE DEGREE TO MEDIAN CONNECTIVITY  
  df_dx_graph <- df_dx_graph %>% 
    left_join(df_dx_sft4 %>% dplyr::select(dx, median_k, r_sq),
              by = "dx") %>% 
    mutate(adj_degree = degree/median_k)
  
# PLOT
  load(paste0("objects/", ctrl_prefix, "_4DX_SFT4_GRAPH.RDS"))
  
  df_dx_graph %>% 
    ggplot(aes(x = adj_degree)) +
    geom_density(aes(fill = dx), alpha = 0.5) +
    facet_wrap( ~ module, scales = "free") +
    ggtitle("Degree distributions by module")
  
  df_dx_graph %>% 
    ggplot(aes(x = hub_score)) +
    geom_density(aes(fill = dx), alpha = 0.5) +
    facet_wrap( ~ module, scales = "free") +
    ggtitle("Hub score distributions by module")
  
# CALCULATE DIFFERENCES USING KRUSKAL-WALLIS & PAIRWISE WILCOX
  
  df_dx_graph_stats <- df_dx_graph %>% 
    group_by(module) %>% 
    nest() %>% 
    
    # KW
    mutate(degree_kw_pval = map(.x = data,
                                .f = ~ kruskal.test(.x$adj_degree, .x$dx)$p.value),
           hub_kw_pval = map(.x = data,
                             .f = ~ kruskal.test(.x$hub_score, .x$dx)$p.value)
    ) %>% 
    unnest(cols = c(degree_kw_pval, hub_kw_pval)) %>% 
    ungroup %>% 
    mutate(degree_kw_qval = p.adjust(degree_kw_pval, method = "bonferroni"),
           hub_kw_qval = p.adjust(hub_kw_pval, method = "bonferroni")) %>% 
    
    # wilcox
    mutate(degree_wilcox_pval = map(.x = data,
                                    .f = ~ pairwise.wilcox.test(.x$adj_degree, .x$dx)$p.value %>% 
                                      as.data.frame %>% 
                                      rownames_to_column("degree_dx1") %>% 
                                      as_tibble() %>% 
                                      pivot_longer(2:ncol(.), names_to = "degree_dx2", values_to = "degree_wilcox_pval") %>% 
                                      filter(!is.na(degree_wilcox_pval))
    ),
    hub_wilcox_pval = map(.x = data,
                          .f = ~ pairwise.wilcox.test(.x$hub_score, .x$dx)$p.value %>% 
                            as.data.frame %>% 
                            rownames_to_column("hub_dx1") %>% 
                            as_tibble() %>% 
                            pivot_longer(2:ncol(.), names_to = "hub_dx2", values_to = "hub_wilcox_pval") %>% 
                            filter(!is.na(hub_wilcox_pval))
    )
    ) %>% 
    unnest(cols = c(degree_wilcox_pval, hub_wilcox_pval)) %>% 
    mutate(degree_wilcox_qval = p.adjust(degree_wilcox_pval, method = "bonferroni"),
           hub_wilcox_qval = p.adjust(hub_wilcox_pval, method = "bonferroni"))
  
  
  df_dx_graph_stats %>% 
    filter(degree_kw_qval < 0.05 & degree_wilcox_qval < 0.05) %>% 
    select(module, degree_dx1, degree_dx2, degree_wilcox_qval) %>% 
    print(n = nrow(.))
  
  


# identify sft where diagnostic groups have similar fit & median k --------

  
# CALCULATE SFT FOR EACH DIAGNOSIS
  powers <- 1:16
  
  df_dx_sft <- df_vsd_regress %>% 
    pivot_wider(id_cols = sample, 
                names_from = ensembl_gene_id, 
                values_from = resids) %>% 
    mutate(dx = gsub(".*_", "", sample), .before = 1) %>% 
    group_by(dx) %>% 
    nest() %>% 
    
    mutate(data_mat = map(.x = data,
                          .f = ~ .x %>%     
                            as.data.frame %>% 
                            column_to_rownames("sample") %>% 
                            as.matrix)) %>% 
    mutate(sft = map(.x = data,
                     .f = ~ pickSoftThreshold(.x, powerVector = powers, verbose = 5)))
  
 
  save(df_dx_sft,
       file = paste0("objects/", ctrl_prefix, "_4DX_SFT.RDS"))
  
  
  
# PLOT FIT METRICS FOR EACH DX  
  
  load(paste0("objects/", ctrl_prefix, "_4DX_SFT.RDS"))
  
  df_dx_sft %>% 
    mutate(sft = map(.x = sft,
                     .f = ~ .x$fitIndices %>% 
                       as_tibble %>% 
                       clean_names %>% 
                       dplyr::select(power, sft_r_sq, mean_k, median_k, max_k) %>% 
                       pivot_longer(2:ncol(.), names_to = "metric", values_to = "value")
                     )
           ) %>% 
    unnest(cols = c(sft)) %>% 
    
    ggplot(aes(x = power, y = value, color = dx)) +
    geom_point(shape = 1, size = 3) +
    facet_wrap(~ metric, scales = "free") +
    ggtitle("At what sft do fit metrics across Dx groups converge?")
  
  
# CALCULATE ADJACENCY AND TOM FOR EACH DIAGNOSTIC GROUP SEPARATELY (SFT = 8)  
  df_dx_wgcna <- df_vsd_regress %>% 
    pivot_wider(id_cols = sample, 
                names_from = ensembl_gene_id, 
                values_from = resids) %>% 
    mutate(dx = gsub(".*_", "", sample), .before = 1) %>% 
    group_by(dx) %>% 
    nest() %>% 
    
    mutate(data_mat = map(.x = data,
                          .f = ~ .x %>%     
                            as.data.frame %>% 
                            column_to_rownames("sample") %>% 
                            as.matrix)) %>% 
    
    mutate(adj = map(.x = data_mat,
                     .f = ~ adjacency(.x, power = 8)
    )) %>% 
    
    mutate(tom = map(.x = adj,
                     .f = ~ TOMsimilarity(.x))
    ) %>% 
    dplyr::select(dx, tom)
  
  save(df_dx_wgcna,
       file = paste0("objects/", ctrl_prefix, "_4DX_SFT8_TOM.RDS"))
  
  
  
# PIVOT TOM INTO IGRAPH FORMAT
  load(paste0("objects/", ctrl_prefix, "_4DX_SFT8_TOM.RDS")) # df_dx_wgcna
  
  df_dx_tom <- df_dx_wgcna %>% 
    
    mutate(tom = map(.x = tom,
                     .f = ~ .x %>% 
                       as.data.frame %>% 
                       as_tibble %>% 
                       dplyr::rename_all( ~ unique(df_vsd_regress$ensembl_gene_id)) %>% 
                       mutate(gene1 = unique(df_vsd_regress$ensembl_gene_id), .before = 1) %>% 
                       pivot_longer(starts_with("ENSG"), names_to = "gene2", values_to = "weight") %>% 
                       
                       left_join(df_modules_filt %>% 
                                   dplyr::rename("gene1" = "ensembl_gene_id",
                                                 "module1" = "module"),
                                 by = "gene1"
                       ) %>% 
                       
                       left_join(df_modules_filt %>% 
                                   dplyr::rename("gene2" = "ensembl_gene_id",
                                                 "module2" = "module"),
                                 by = "gene2"
                       ) %>% 
                       
                       filter(module1 == module2) %>% 
                       dplyr::select(-module2) %>% 
                       dplyr::rename("module" = "module1")
    )
    
    )
  
  save(df_dx_tom,
       file = paste0("objects/", ctrl_prefix, "_4DX_SFT8_TOM.RDS")) 
  
  
# GENERATE GRAPH FOR EACH MODULE IN EACH DX GROUP  
  load(paste0("objects/", ctrl_prefix, "_4DX_SFT8_TOM.RDS")) # df_dx_tom
  
  df_dx_graph <- df_dx_tom %>% 
    unnest(cols = c(tom)) %>% 
    group_by(dx, module) %>% 
    nest() %>% 
    mutate(graph = map(.x = data,
                       .f = ~ graph_from_data_frame(.x, directed = FALSE))) %>% 
    dplyr::select(dx, module, graph) %>% 
    mutate(dx = factor(dx, levels = c("control", "bipolar", "mdd", "schizo"))) %>% 
    arrange(dx, module)
  
  save(df_dx_graph,
       file = paste0("objects/", ctrl_prefix, "_4DX_SFT8_GRAPH.RDS"))
  
# CALCULATE DEGREE AND HUB SCORE
  load(paste0("objects/", ctrl_prefix, "_4DX_SFT8_GRAPH.RDS")) # df_dx_graph
  
  df_dx_graph <- df_dx_graph %>% 
    mutate(degree = map(.x = graph,
                        .f = ~ graph.strength(.x)
    ),
    hub_score = map(.x = graph,
                    .f = ~ hub_score(.x)$vector
    )
    )  %>% 
    unnest(cols = c(degree, hub_score))
  
  save(df_dx_graph,
       file = paste0("objects/", ctrl_prefix, "_4DX_SFT8_GRAPH.RDS"))
  

# PLOT
  df_dx_graph %>% 
    ggplot(aes(x = degree)) +
    geom_density(aes(fill = dx), alpha = 0.5) +
    facet_wrap( ~ module, scales = "free") +
    ggtitle("Degree distributions by module")
  
  df_dx_graph %>% 
    ggplot(aes(x = hub_score)) +
    geom_density(aes(fill = dx), alpha = 0.5) +
    facet_wrap( ~ module, scales = "free") +
    ggtitle("Hub score distributions by module")
  
# CALCULATE DIFFERENCES USING KRUSKAL-WALLIS & PAIRWISE WILCOX
  
  df_dx_graph_stats <- df_dx_graph %>% 
    group_by(module) %>% 
    nest() %>% 
    
    # KW
    mutate(degree_kw_pval = map(.x = data,
                                .f = ~ kruskal.test(.x$degree, .x$dx)$p.value),
           hub_kw_pval = map(.x = data,
                             .f = ~ kruskal.test(.x$hub_score, .x$dx)$p.value)
    ) %>% 
    unnest(cols = c(degree_kw_pval, hub_kw_pval)) %>% 
    ungroup %>% 
    mutate(degree_kw_qval = p.adjust(degree_kw_pval, method = "bonferroni"),
           hub_kw_qval = p.adjust(hub_kw_pval, method = "bonferroni")) %>% 
    
    # wilcox
    mutate(degree_wilcox_pval = map(.x = data,
                                    .f = ~ pairwise.wilcox.test(.x$degree, .x$dx)$p.value %>% 
                                      as.data.frame %>% 
                                      rownames_to_column("degree_dx1") %>% 
                                      as_tibble() %>% 
                                      pivot_longer(2:ncol(.), names_to = "degree_dx2", values_to = "degree_wilcox_pval") %>% 
                                      filter(!is.na(degree_wilcox_pval))
    ),
    hub_wilcox_pval = map(.x = data,
                          .f = ~ pairwise.wilcox.test(.x$hub_score, .x$dx)$p.value %>% 
                            as.data.frame %>% 
                            rownames_to_column("hub_dx1") %>% 
                            as_tibble() %>% 
                            pivot_longer(2:ncol(.), names_to = "hub_dx2", values_to = "hub_wilcox_pval") %>% 
                            filter(!is.na(hub_wilcox_pval))
    )
    ) %>% 
    unnest(cols = c(degree_wilcox_pval, hub_wilcox_pval)) %>% 
    mutate(degree_wilcox_qval = p.adjust(degree_wilcox_pval, method = "bonferroni"),
           hub_wilcox_qval = p.adjust(hub_wilcox_pval, method = "bonferroni"))
  
  
  df_dx_graph_stats %>% 
    filter(degree_kw_qval < 0.05 & degree_wilcox_qval < 0.05) %>% 
    dplyr::select(module, degree_dx1, degree_dx2, degree_wilcox_qval) %>% 
    print(n = nrow(.))
  


   
  
# Louvain community detection ---------------------------------------------


# RUN LOUVAIN CLUSTERING OF SAMPLES BASED ON GENES IN EACH MODULE
  df_louvain <- df_vsd_regress %>% 
    left_join(df_modules_filt) %>% 
    dplyr::select(sample, module, ensembl_gene_id, resids) %>% 
    group_by(module) %>% 
    nest() %>% 
    arrange(module) %>% 
    
    mutate(data = map(.x = data,
                      .f = ~ .x %>% 
                        pivot_wider(id_cols = ensembl_gene_id, names_from = sample, values_from = resids) %>% 
                        correlate() %>% 
                        dplyr::rename("sample1" = "term") %>% 
                        pivot_longer(2:ncol(.), names_to = "sample2", values_to = "weight") %>% 
                        filter(!is.na(weight)) %>% 
                        mutate(weight = abs(weight))
                      )
           ) %>% 
    
    mutate(graph = map(.x = data,
                       .f = ~ graph_from_data_frame(.x, directed = FALSE)
                       )
           ) %>% 
    
    mutate(louvain = map(.x = graph,
                         .f = ~ cluster_louvain(.x)))
  
  save(df_louvain,
       file = paste0("objects/", ctrl_prefix, "_4DX_SFTK_LOUVAIN.RDS"))
  
# DETERMINE MODULE MEMBERSHIP FOR EACH SAMPLE 
  df_louvain_res <- df_louvain %>% 
    dplyr::select(module, louvain) %>% 
    mutate(louvain = map(.x = louvain,
                         .f = ~ tibble(
                           
                           sample = .x$names,
                           cluster = .x$membership
                           
                         )
    )
    ) %>% 
    unnest(cols = c(louvain))
  
# PLOT
  df_louvain_res %>% 
    mutate(dx = gsub(".*_", "", sample) %>% factor(levels = c("control", "bipolar", "mdd", "schizo"))) %>% 
    group_by(module, dx) %>% 
    dplyr::count(cluster) %>% 
    dplyr::summarise(cluster = as.factor(cluster), 
                     perc = 100*n/sum(n)) %>% 
    
    ggplot(aes(x = cluster, y = perc, fill = dx)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~module) +
    geom_hline(aes(yintercept = 50), lty = 2, color = "black") +
    ylim(c(0, 100)) +
    labs(y = "percent samples in group") +
    ggtitle("Louvain")
  
    
# PERMS
  
  df_louvain_res %>% 
    group_by(module) %>% 
    nest() %>% 
    
    mutate(perm_clust = map(.x = data,
                            .f = ~ sample(
                              x = .x$cluster %>% unique,
                              replace = TRUE,
                              size = length(.x$cluster),
                              prob = .x %>%
                                dplyr::count(cluster) %>%
                                summarise(prop = n/sum(n)) %>%
                                pull(prop)
                            )
    )
    ) %>% 
    unnest(cols = c(data, perm_clust))
  
  props <- df_louvain_res %>% 
    group_by(module) %>% 
    nest() %>% 
    
    pull(data) %>%  .[[2]] %>%
    
    dplyr::count(cluster) %>%
    summarise(prop = n/sum(n)) %>% pull(prop)
  
  cluster <- df_louvain_res %>% 
    group_by(module) %>% 
    nest() %>% 
    
    pull(data) %>%  .[[2]] %>% pull(cluster) %>% unique
    
    sample(cluster, size = length(cluster), replace = TRUE, prob = props)
  
    df_louvain_res %>% 
      group_by(module) %>% 
      nest() %>% 
      
      pull(data) %>%  .[[2]] %>% 
      
      mutate(cluster = factor(cluster)) %>% 
      mutate(rand = sample(cluster, n(), replace = TRUE, prob = c(0.45, 0.45, 0.05, 0.05)))
      
  df_louvain_perms <- tibble()
  for (i in 0:max(as.numeric(as.character(df_louvain_res$module)))) {
    
    df_mod <- df_louvain_res %>% 
      filter(module == 1) %>% 
      mutate(dx = gsub(".*_", "", sample))
    
    df_tmp <- 1000 %>% rerun(
      
      sample(df_mod %>% pull(cluster) %>% unique,
             size = nrow(df_mod),
             replace = TRUE,
             prob = df_mod %>% 
               dplyr::count(cluster) %>% 
               summarise(prop = n/sum(n)) %>% 
               pull(prop))      
      
    ) %>% 
      bind_cols %>% 
      rename_all(~paste0("perm_", str_remove(.x, "..."))) %>% 
      bind_cols(df_mod, .)
    
    df_louvain_perms <- df_louvain_perms %>% bind_rows(df_tmp)
    
  }
  
  df_louvain_perms_pvals <- df_louvain_perms %>% 
    pivot_longer(contains("perm"), names_to = "perm_no", values_to = "perm_group") %>% 
    dplyr::count(module, dx, perm_no, perm_group) %>% 
    dplyr::rename("perm_n" = "n") %>% 
    group_by(module, dx) %>%
    left_join(df_louvain_perms %>% 
                dplyr::count(module, dx, cluster) %>% 
                dplyr::rename("actual_n" = "n")) %>% 
    filter(perm_group == cluster) %>% 
    mutate(extreme = ifelse(perm_n > actual_n, 1, 0)) %>% 
    dplyr::count(cluster, extreme) %>% 
    filter(extreme == 1) %>% 
    dplyr::select(-extreme) %>% 
    ungroup %>% 
    mutate(perm_p_val = n/1000,
           perm_q_val = p.adjust(perm_p_val, method = "BH"))
  
  df_louvain_perms_pvals %>% 
    filter(perm_p_val < 0.05)
  
  
  

