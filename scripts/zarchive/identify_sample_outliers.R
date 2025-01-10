


# libraries ---------------------------------------------------------------


    library(tidyverse)
    library(janitor)
    library(plotly)
    library(RColorBrewer)
    library(ggdendro)
    library(microplot)
    library(tidymodels)            
    library(purrr)
    library(grid)
    library(dendextend)
    library(themis)
    library(tidytext)


# set theme for plots -----------------------------------------------------

  theme_set(theme_bw() +
              theme(plot.title = element_text(size = 18),
                    axis.title = element_text(size = 15),
                    axis.text = element_text(size = 12),
                    strip.text = element_text(size = 15),
                    legend.title = element_text(size = 15),
                    legend.text = element_text(size = 12)))


# data --------------------------------------------------------------------


# NORMALIZED AND FILTERED COUNTS
    # 14k genes, vsd regressed, RIN_GC_race_libraryBatch regressed
    df_vsd_resids <- read.csv("outputs/20220613_185samples_14kgenes_vsd_RIN_mappedPerc_GC_race_resids.csv") %>% 
      as_tibble() %>% 
      dplyr::select(-X) %>% 
      rename_at(2:ncol(.), ~str_remove(.x, "X"))
    
    # transpose
    df_vsd_resids_t <- df_vsd_resids %>%
      
      # transpose & convert to tibble
      dplyr::select(-1) %>%
      t %>%
      as.data.frame %>%
      as_tibble %>%
      
      # rename columns with gene id from original df
      rename_all(~df_vsd_resids %>% pull(ensembl_gene_id)) %>%
      
      # add sample id col using colnames of original df (to maintain order)
      dplyr::mutate(sample = df_vsd_resids %>% dplyr::select(-1) %>% colnames, .before = 1)

# COVARIATES
    load("data/covariates/185_all_covariates_clean.Rdata") # df_covariates_all_clean, df_covariates_all_clean_numeric

    

# gene tree ---------------------------------------------------------------


    # cluster norm counts
    sample_tree <- df_vsd_resids_t %>% 
      as.data.frame %>%
      column_to_rownames("sample") %>% 
      dist %>% 
      hclust
    
    # plot
    plot(sample_tree, cex = 0.6)
    

# PCA ---------------------------------------------------------------------

    
    # SAMPLE PCA
    
    pca <- prcomp(df_vsd_resids %>% 
                    as.data.frame %>%
                    column_to_rownames("ensembl_gene_id"))
    
    df_pca <- pca$rotation %>% 
      as.data.frame %>% 
      rownames_to_column("sample") %>% 
      as_tibble() %>% 
      left_join(df_covariates_all_clean) %>% 
      mutate(library_batch = as.factor(library_batch)) %>% 
      clean_names()
    
    df_var <- summary(pca) %>% 
      .$importance %>% 
      as.data.frame %>% 
      rownames_to_column("metric") %>% 
      as_tibble() %>% 
      pivot_longer(2:ncol(.), values_to = "value", names_to = "PC") %>% 
      mutate(PC = str_remove(PC, "PC") %>% as.numeric()) %>% 
      pivot_wider(id_cols = PC, names_from = metric, values_from = value) %>% 
      clean_names()
    
    p_pca <- df_pca %>% 
      ggplot(aes(x = pc1, y = pc2)) +
      geom_point(color = "black", size = 2.75) +
      geom_point(aes(color = new_dx, text = sample), size = 2) +
      
      # label samples
      # geom_text_repel(aes(label = sample)) +
      
      # add origin lines
      geom_vline(xintercept = 0, color = "black", lty = 2) +
      geom_hline(yintercept = 0, color = "black", lty = 2) +
      
      # give % variance in axis labels
      xlab(paste0("PC1: ", 
                  df_var %>% filter(pc == 1) %>% pull(proportion_of_variance) *100 %>% round(2),
                  "% of variance")) +
      ylab(paste0("PC2: ", 
                  df_var %>% filter(pc == 2) %>% pull(proportion_of_variance) *100 %>% round(2),
                  "% of variance"))
    
    ggplotly(p_pca, tooltip = "text")
    
    # plot eigenvalues
    df_var %>% # filter(pc != 1) %>% 
      ggplot() +
      geom_col(aes(x = pc, y = proportion_of_variance*100), fill = "lightblue", color = "black") +
      
      # error bars
      # geom_segment(aes(x = pc, y = proportion_of_variance*100 - standard_deviation,
      #                  xend = pc, yend = proportion_of_variance*100 + standard_deviation)) +
      
      labs(y = "percent variance") +
      scale_x_continuous(breaks = seq(5, max(df_var$pc), 5)) +
      theme_bw()



# covariate-PC correlations -----------------------------------------------


    # RUN CORREALTIONS
    pcs <- df_pca %>% dplyr::select(contains("pc")) %>% colnames
    covariates <- df_covariates_all_clean_numeric %>% dplyr::select(2:ncol(.)) %>% colnames
    df_pc_covar_cor <- tibble()
    
    for (pc in pcs) {
      
      for (var in covariates) {
        
        cor <- cor.test(df_pca[[pc]], df_covariates_all_clean_numeric[[var]])
        estimate <- cor$estimate
        p <- cor$p.value
        
        df_tmp <- tibble(pc = pc,
                         covariate = var,
                         cor = estimate,
                         p_val = p)
        
        df_pc_covar_cor <- df_pc_covar_cor %>% bind_rows(df_tmp)
        
      }
      
    }
    
    df_pc_covar_cor <- df_pc_covar_cor %>% 
      mutate(pc = str_remove(pc, "pc") %>% as.numeric()) %>% 
      left_join(df_var %>% 
                  pivot_longer(2:ncol(.), names_to = "metric", values_to = "value") %>% 
                  filter(metric == "proportion_of_variance"), by = "pc") %>% 
      dplyr::select(-metric) %>% 
      relocate(value, .before = covariate) %>% 
      dplyr::rename("pc_perc_var" = "value") %>% 
      mutate(pc_perc_var = 100*pc_perc_var) %>% 
      mutate(p_adj = p.adjust(p_val, method = "BH")) # %>% write.csv("outputs/pc_14k_genes_vsd_regress_RIN_race_GC_libraryBatch_significant_covariate_correlations.csv")
    
    df_pc_covar_cor %>% filter(p_adj < 0.05 & pc_perc_var > 0)
    
    # PLOT HEATMAP
    df_pc_covar_cor %>% 
      filter(pc_perc_var > 1) %>% 
      mutate(p_adj = ifelse(p_adj == 1, NA, p_adj)) %>% 
      mutate(labels = ifelse(p_adj < 0.05, formatC(p_adj, format = "e", digits = 2), "")) %>% 
      
      ggplot(aes(x = factor(pc), y = covariate)) +
      geom_tile(aes(fill = -log10(p_adj))) +
      geom_text(aes(label = labels), color = "white") +
      scale_fill_gradientn(colors = brewer.pal(9, "Reds"), na.value = "white") +
      labs(x = "PC (> 1% explained var)") +
      theme(axis.text.y = element_text(size = 12))
    
    
    
    
# z-score PCA -------------------------------------------------------------

    
    # FIND Z-SCORE AT SIGNIFICANT PCS
    sig_pcs <- df_var %>% filter(proportion_of_variance > 0.01) %>% pull(pc) %>% as.character()
    
    df_pc_zscore <- df_pca %>% 
      dplyr::select(sample, starts_with("pc")) %>% 
      pivot_longer(2:ncol(.), values_to = "value", names_to = "pc") %>% 
      group_by(pc) %>% 
      mutate(z_score = (value - mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE)) %>% 
      mutate(pc = str_remove_all(pc, "pc")) %>% 
      filter(pc %in% sig_pcs)
      
    # DENDROGRAM TO CLUSTER
    pc_zscore_tree <- df_pc_zscore %>% 
      pivot_wider(id_cols = sample, names_from = pc, values_from = z_score) %>% 
      as.data.frame %>%
      column_to_rownames("sample") %>% 
      dist %>% 
      hclust
    
    pc_zscore_dendro <- as.dendrogram(pc_zscore_tree)
    labels_pc_zscore <- labels(pc_zscore_dendro)
    
    # PLOT HEATMAP
    df_pc_zscore %>% 
      
      ggplot(aes(x = pc, y = sample)) +
      geom_tile(aes(fill = z_score)) +
      scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), na.value = "white") +
      scale_x_discrete(limits = sig_pcs) +
      scale_y_discrete(limits = labels_pc_zscore) +
      # labs(x = "", y = "") +
      theme(axis.text.y = element_text(size = 8))
    
    # SEE IF ANY SAMPLES POP UP AS OUTLIERS ON A LOT OF PCS
    df_pc_zscore %>% 
      ungroup %>% 
      filter(abs(z_score) > 3) %>% 
      dplyr::count(sample) %>% 
      arrange(pc) %>% print(n = nrow(.))
    
    
    

# correlation matrix ------------------------------------------------------


        m_sample_cor <- cor(df_vsd_resids[2:ncol(df_vsd_resids)], method = "pearson")
        diag(m_sample_cor) <- NA
        
        df_sample_cor <- m_sample_cor %>% 
          as.data.frame %>% 
          rownames_to_column("sample1") %>% 
          as_tibble %>% 
          pivot_longer(2:ncol(.), names_to = "sample2", values_to = "pearsons_r")
        
        df_sample_cor_zscore <- df_sample_cor %>% 
          mutate(z_score = (pearsons_r - mean(pearsons_r, na.rm = TRUE))/sd(pearsons_r, na.rm = TRUE))
    
        m_sample_cor_zscore <- df_sample_cor_zscore %>% 
          pivot_wider(id_cols = sample1, names_from = sample2, values_from = z_score) %>% 
          as.data.frame %>% 
          column_to_rownames("sample1") %>% 
          as.matrix

    # GEOM TILE & GGDENDRO HEATMAP
        
        sample_dendro <- as.dendrogram(sample_tree)
        labels <- labels(sample_dendro)
        
        # plot heatmap
        p_hm <- df_sample_cor_zscore %>% 
          
          ggplot(aes(x = sample1, y = sample2)) +
          geom_tile(aes(fill = z_score)) +
          scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), na.value = "white") +
          scale_x_discrete(limits = labels) +
          scale_y_discrete(limits = labels) + 
          labs(x = "", y = "") +
          theme(# axis.text.x = element_text(angle = 90, hjust = 0.9),
                axis.text = element_blank(),
                # axis.text = element_text(size = 8),
                axis.ticks = element_blank())
        
        # plot dendrogram
        ddata <- dendro_data(sample_dendro, type = "rectangle")
        ddata$segments["yend"][ddata$segments["yend"] == 0] <- 30
        
        p_dendro <- ggplot(segment(ddata)) + 
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
          theme_collapse() +
          theme(panel.border = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, "pt"),
                plot.margin = margin(0, 0, 0, 0, "pt"))

        # annotation bar
        library_batch <- df_covariates_all_clean %>% 
          mutate(library_batch = as.factor(library_batch)) %>% 
          ggplot(aes(x = sample, y = 0.1)) +
          geom_point(aes(color = library_batch), shape = 15, size = 1.3, show.legend = F) +
          scale_x_discrete(limits = labels) +
          theme_collapse()
        
        dx <- df_covariates_all_clean %>% 
          ggplot(aes(x = sample, y = 0.1)) +
          geom_point(aes(color = new_dx), shape = 15, size = 1.3, show.legend = F) +
          scale_x_discrete(limits = labels) +
          theme_collapse()
        manner <- df_covariates_all_clean %>% 
          ggplot(aes(x = sample, y = 0.1)) +
          geom_point(aes(color = manner), shape = 15, size = 1.3, show.legend = F) +
          scale_x_discrete(limits = labels) +
          scale_color_manual(values = c("grey", "green", "red", "blue")) +
          theme_collapse()
        opioids_sedatives <- df_covariates_all_clean %>% 
          ggplot(aes(x = sample, y = 0.1)) +
          geom_point(aes(color = opioids_sedatives), shape = 15, size = 1.3, show.legend = F) +
          scale_x_discrete(limits = labels) +
          theme_collapse()
        
        
        # put together
        grid.newpage()
        print(p_hm,
              vp = viewport(x = 0.5, y = 0.47, width = 0.8, height = 0.8))
        print(p_dendro,
              vp = viewport(x = 0.468, y = 0.92, width = 0.74, height = 0.13))
        
        print(library_batch,
              vp = viewport(x = 0.465, y = 0.1, width = 0.6795, height = 0.3))
        print(dx,
              vp = viewport(x = 0.465, y = 0.1, width = 0.6795, height = 0.3))
        print(manner,
              vp = viewport(x = 0.465, y = 0.09, width = 0.6795, height = 0.3))
        print(opioids_sedatives,
              vp = viewport(x = 0.465, y = 0.08, width = 0.6795, height = 0.3))
        


        
# sample clusters ---------------------------------------------------------


        # PLOT HEATMAP WITH CLUSTERED DENDRO
        
            p_hm <- df_sample_cor_zscore %>%

              ggplot(aes(x = sample1, y = sample2)) +
              geom_tile(aes(fill = z_score)) +
              scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), na.value = "white") +
              scale_x_discrete(limits = labels_nobd) +
              scale_y_discrete(limits = labels_nobd) +
              labs(x = "", y = "") +
              theme(axis.text.x = element_text(angle = 90, hjust = 0.9),
                    # axis.text = element_blank(),
                    axis.text = element_text(size = 8),
                    axis.ticks = element_blank())

            # plot dendrogram
            sample_dendro <- as.dendrogram(sample_tree) %>% set("branches_k_color", k = 3)
            labels <- labels(sample_dendro)
            ddata <- dendro_data(sample_dendro, type = "rectangle")
            df_clust <- cutree(sample_tree, k = 3) %>% # 3 k-means clusters from hclust data
              enframe("label", "cluster") %>% 
              mutate(cluster = as.factor(cluster))
            
            ddata[["labels"]] <- merge(ddata[["labels"]], df_clust, by = "label")
            ddata$segments["yend"][ddata$segments["yend"] == 0] <- 30

            p_dendro <- ggplot() +
              geom_segment(data = segment(ddata), 
                           aes(x = x, y = y, xend = xend, yend = yend))+
              geom_segment(data = ddata$segments %>%
                             filter(yend == 30) %>%
                             left_join(ddata$labels, by = "x"), 
                           aes(x = x, y = y.x, xend = xend, yend = yend, color = cluster)) +
              theme_collapse() +
              theme(panel.border = element_blank(),
                    axis.ticks = element_blank(),
                    axis.ticks.length = unit(0, "pt"),
                    plot.margin = margin(0, 0, 0, 0, "pt"),
                    legend.position = "none")

            # put together
            grid.newpage()
            print(p_hm,
                  vp = viewport(x = 0.5, y = 0.47, width = 0.8, height = 0.8))
            print(p_dendro,
                  vp = viewport(x = 0.468, y = 0.92, width = 0.74, height = 0.13))


            

# explore covariates driving clusters -------------------------------------



      # EXPLORE COVARIATES IN EACH CLUSTER
            
            # separate numeric and categorical variables
            numeric_cov <- df_covariates_all_clean %>% dplyr::select_if(is.numeric) %>% colnames
            numeric_cov <- numeric_cov[! numeric_cov %in% c("multi_nomial_dx", "library_batch", "rn_aextraction_batch")]
            categorical_cov <- c(df_covariates_all_clean %>% dplyr::select_if(is.character) %>% colnames,
                                 "library_batch", "rn_aextraction_batch")
            
            # plots

                # numeric (histogram)
                df_numeric_cov_clust_plot <- df_clust %>% 
                  dplyr::rename("sample" = "label") %>% 
                  left_join(df_clust %>% dplyr::count(cluster)) %>% 
                  right_join(df_covariates_all_clean %>%
                               dplyr::select(sample, contains(numeric_cov)) %>% 
                               pivot_longer(2:ncol(.), 
                                            values_to = "value", 
                                            names_to = "covariate")) %>%
                  filter(!(covariate %in% c("hallucinogens", "cocaine")))
                
                # histogram
                df_numeric_cov_clust_plot %>% 
    
                  ggplot(aes(x = value)) +
                  geom_density(aes(color = cluster, fill = cluster), alpha = 0.3, position = "dodge") +
                  facet_wrap(~covariate, scales = "free") +
                  labs(x = "") +
                  theme(strip.text = element_text(size = 12))
            
              # categorial (frequency)
                df_categorical_cov_clust_plot <- df_clust %>% 
                  dplyr::rename("sample" = "label") %>% 
                  left_join(df_clust %>% dplyr::count(cluster)) %>% 
                  right_join(df_covariates_all_clean %>%
                               dplyr::select(sample, contains(categorical_cov)) %>% 
                               mutate_all(~ifelse(.x == "", "NA", .x)) %>% 
                               mutate_all(~as.factor(as.character(.x))) %>% 
                               pivot_longer(2:ncol(.), 
                                            values_to = "value", 
                                            names_to = "covariate")) %>%
                  filter(!(covariate %in% c("hallucinogens", "cocaine")))
                
                # bar plot
                df_categorical_cov_clust_plot %>% 
                  dplyr::select(-value) %>% 
                  left_join( df_categorical_cov_clust_plot %>% 
                               dplyr::count(cluster, n, covariate, value), 
                             by = c("cluster", "n", "covariate")) %>% 
                  mutate(percent = 100 * nn/n) %>% 
                  mutate(value = as.factor(as.character(value))) %>% 
                  
                  ggplot(aes(x = percent, y = str_wrap(value, width = 15))) +
                  geom_bar(aes(color = cluster, fill = cluster), alpha = 0.3, position = "dodge", stat = "identity") +
                  facet_wrap(~covariate, scales = "free") +
                  labs(y = "") +
                  theme(strip.text = element_text(size = 10))
                
                # dot plot
                df_clust_num <- df_clust_covariate_tests %>% 
                  dplyr::count(cluster) %>% 
                  dplyr::rename("cluster_n" = "n") %>% 
                  mutate(cluster = as.factor(cluster))
                
                df_categorical_cov_clust_plot %>% 
                  dplyr::count(cluster, covariate, value) %>% 
                  left_join(df_clust_num) %>% 
                  mutate(percent = 100 * n/cluster_n) %>% 
                  
                  ggplot(aes(x = cluster, y = str_wrap(value, width = 15))) +
                  geom_point(aes(size = percent + 1), color = "black") +
                  geom_point(aes(color = percent, size = percent)) +
                  facet_wrap(~covariate, scales = "free") +
                  scale_color_gradientn(colors = brewer.pal(9, "Reds")[3:9]) +
                  labs(y = "") +
                  theme(strip.text = element_text(size = 10))
                
                
            
          # STATISTICAL TESTS

            df_clust_covariate_tests <- df_clust %>%
              dplyr::rename("sample" = "label") %>% 
              right_join(df_covariates_all_clean_numeric)  %>%
              mutate(cluster = as.numeric(cluster))

            # cluster ~ covariate
            df_clust_covariate_res <- numeric_cov %>% map_df(~tibble(

              covariate = .x,
              p_val = aov(df_clust_covariate_tests[["cluster"]] ~ df_clust_covariate_tests[[.x]]) %>%
                tidy() %>%
                pull(p.value) %>%
                .[1]

            ))
            
            df_clust_covariate_res %>% filter(p_val < 0.1)
            

            # covariate ~ cluster
            df_covariate_cluster_res <- numeric_cov %>% map_df(~tibble(
              
              covariate = .x,
              p_val = aov(df_clust_covariate_tests[[.x]] ~ df_clust_covariate_tests[["cluster"]]) %>%
                tidy() %>%
                pull(p.value) %>%
                .[1]
              
            ))
            
            df_covariate_cluster_res %>% filter(p_val < 0.1)
            
            # Chi-sq test
            df_categorical_cov_clust_chisq <- df_categorical_cov_clust_plot %>% 
              dplyr::count(cluster, covariate, value) %>% 
              group_by(covariate) %>% 
              nest() %>% 
              mutate(data = map(data, 
                                ~.x %>%
                                  
                                  pivot_wider(
                                    id_cols = value,
                                    names_from = cluster,
                                    values_from = n
                                  ) %>% 
                                  
                                  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>% 
                                  as.data.frame() %>% 
                                  column_to_rownames("value") %>% 
                                  chisq.test() %>% 
                                  .$p.value)) %>% 
              
              unnest(cols = c(data)) %>%
              ungroup() %>% 
              dplyr::rename("p_val" = "data") %>% 
              mutate(p_adj = p.adjust(p_val, method = "BH"))

            df_categorical_cov_clust_chisq %>% filter(p_val < 0.05)

            

# explore genes driving clusters ------------------------------------------

            
            
      # join cluster and expression data
            df_clust_vsd_resids <- df_clust %>% 
              dplyr::rename("sample" = "label") %>% 
              left_join(df_vsd_resids_t) %>% 
              dplyr::select(-sample)
            
      # HISTOGRAM           
            df_clust_vsd_resids %>% 
              pivot_longer(2:ncol(.), names_to = "ensembl_gene_id", values_to = "residual") %>% 
              
              ggplot(aes(x = log10(residual + 1))) +
              geom_density(aes(fill = cluster, y = stat(density)), 
                           position = "dodge", alpha = 0.5) +
              scale_y_continuous(labels = scales::percent_format())
            
      # MODEL
            
            # split data
            set.seed(20220615)
            clust_split <- df_clust_vsd_resids %>% initial_split(strata = cluster)
            
            clust_train <- training(clust_split)
            clust_test <- testing(clust_split) 
            
            clust_folds <- vfold_cv(clust_train, strata = cluster, v = 100)
            
            clust_train %>% dplyr::count(cluster)
            
            # define model
            svm_spec <- svm_linear(mode = "classification")
            
            # define recipe
            clust_rec <- recipe(cluster ~ ., data = clust_train) %>% 
              step_nzv() %>% 
              step_smote() %>% 
              step_normalize()
            
            # create workflow (combining recipte & model spec)
            clust_wf <- workflow(clust_rec, svm_spec)
            
            # run model on 100 resamples
            doParallel::registerDoParallel()
            clust_metrics <- metric_set(accuracy, sens, spec)
            
            clust_rs <- fit_resamples(clust_wf,
                                      resamples = clust_folds,
                                      metrics = clust_metrics)
            collect_metrics(clust_rs)
            
            # fun on test data
            final_rs <- last_fit(clust_wf,
                                 clust_split,
                                 metrics = clust_metrics)
            collect_metrics(final_rs)
            
            # plot confusion matrix
            collect_predictions(final_rs) %>%
              count(cluster, .pred_class) %>% 
              
              ggplot(aes(x = cluster, y = .pred_class)) +
              geom_tile(aes(fill = cluster)) +
              geom_text(aes(label = n), size = 12) +
              labs(y = "prediction", x = "truth") +
              theme(axis.text = element_text(size = 15))
            
            # extract final fitted model
            final_fitted <- extract_workflow(final_rs)
            
            # select driving genes
            load("objects/ensembl_id_to_gene_symbol.Rdata")
            tidy(final_fitted) %>% 
              dplyr::rename("ensembl_gene_id" = "term", "coefficient" = "estimate") %>% 
              left_join(df_ensembl_to_symbol) %>% 
              mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
              slice_max(abs(coefficient), n = 30) %>% 
              
              ggplot(aes(x = abs(coefficient), y = reorder(gene_symbol, abs(coefficient)), fill = coefficient > 0)) +
              geom_col() +
              labs(y = "gene")
            
     # FIND GENES THAT HAVE HIGHER/LOWER OVERALL EXPRESSION IN ONE CLUSTER VS ANOTHER
            
            df_clust_zscore <- df_clust_vsd_resids %>% 
              pivot_longer(2:ncol(.), values_to = "value", names_to = "ensembl_gene_id") %>% 
              group_by(ensembl_gene_id) %>% 
              mutate(z_score = (value - mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE))
            
            df_clust_zscore_outliers <- df_clust_zscore %>% 
              filter(z_score > 3) %>% 
              left_join(df_ensembl_to_symbol) %>% 
              mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
              ungroup()
            
            clust1_outliers <- df_clust_zscore_outliers %>% filter(cluster == 1) %>% pull(ensembl_gene_id)
            clust2_outliers <- df_clust_zscore_outliers %>% filter(cluster == 2) %>% pull(ensembl_gene_id)
            clust3_outliers <- df_clust_zscore_outliers %>% filter(cluster == 3) %>% pull(ensembl_gene_id)
            
            # venn diagram
            library(BioVenn)
            par(mar=c(1, 1, 1, 1))
            draw.venn(clust1_outliers, clust2_outliers, clust3_outliers,
                      xtitle = "cluster 1", ytitle = "cluster 2", ztitle = "cluster 3",
                      x_c = "#F8766D", y_c = "#00BA38", z_c = "#619CFF",
                      nr_c = "black")
            
            
            # functional enrichment on outliers
            
                source("functions/f_top_GO_modules.R")
            
                # pull unique gene
                clust1_outliers_only <- clust1_outliers[! clust1_outliers %in% c(clust2_outliers, clust3_outliers)]
                clust2_outliers_only <- clust2_outliers[! clust2_outliers %in% c(clust1_outliers, clust3_outliers)]
                clust3_outliers_only <- clust3_outliers[! clust3_outliers %in% c(clust1_outliers, clust2_outliers)]
                
                # get gene universe
                l_gene_universe <- df_vsd_resids$ensembl_gene_id
                
                # get gene info from ensembl database
                hsapiens_ensembl <- useEnsembl(biomart = "genes",
                                               dataset = "hsapiens_gene_ensembl",
                                               version = 88) # GRCh38.87
                
                # map each gene to GO term
                m_go_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "namespace_1003"),
                                  filters = "ensembl_gene_id",
                                  values = l_gene_universe,
                                  mart = hsapiens_ensembl)
                
                # build annotation list
                l_gene_2_GO <- m_go_ids[,c(1,3)] %>% unstack
                
                # cluster 1
                clust1_outlier_go <- f_top_GO_modules(l_gene_module = clust1_outliers_only,
                                 l_gene_universe = l_gene_universe,
                                 m_go_ids = m_go_ids,
                                 l_gene_2_GO = l_gene_2_GO,
                                 go_ontology = "BP")
                
                # cluster 2
                clust2_outlier_go <- f_top_GO_modules(l_gene_module = clust2_outliers_only,
                                                      l_gene_universe = l_gene_universe,
                                                      m_go_ids = m_go_ids,
                                                      l_gene_2_GO = l_gene_2_GO,
                                                      go_ontology = "BP")
                
                # clsuter 3
                clust3_outlier_go <- f_top_GO_modules(l_gene_module = clust3_outliers_only,
                                                      l_gene_universe = l_gene_universe,
                                                      m_go_ids = m_go_ids,
                                                      l_gene_2_GO = l_gene_2_GO,
                                                      go_ontology = "BP")
                
                # combine results
                load("objects/df_go_full_definitions.Rdata")
                clust_outlier_go <- clust1_outlier_go %>% 
                  mutate(cluster = 1, .before = 1) %>% 
                  bind_rows(clust2_outlier_go %>% 
                              mutate(cluster = 2, .before = 1)) %>% 
                  bind_rows(clust3_outlier_go %>% 
                              mutate(cluster = 3, .before = 1)) %>% 
                  clean_names() %>% 
                  dplyr::select(-term) %>% 
                  left_join(df_go_defs) %>%
                  dplyr::select(-definition) %>%
                  unique
                
                # plot
                clust_outlier_go %>% 
                  group_by(cluster) %>% 
                  arrange(cluster, p_adj) %>% 
                  top_n(n = 10) %>% 
                  
                  ggplot(aes(x = -log10(p_adj), y = reorder_within(str_wrap(term, width = 15), -log10(p_adj), cluster))) +
                  geom_point(aes(color = -log10(p_adj), size = significant/annotated)) +
                  facet_wrap(~cluster, scales = "free_y") +
                  geom_vline(aes(xintercept = -log10(0.05)), lty = 2) +
                  scale_y_discrete(labels = function(x) (gsub("*_.", "", x))) +
                  labs(y = "pathway")
                
                
                
                
                
                
                
                
                
            
            
            
    