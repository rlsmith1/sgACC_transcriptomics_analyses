



# libraries ---------------------------------------------------------------


    library(tidyverse)
    library(purrr)
    library(janitor)
    library(naniar)
    library(Hmisc)
    library(DESeq2)



# data --------------------------------------------------------------------



    # COUNTS
    df_raw_counts <- read.table("data/185samples_allGenes_auto-PARs_noSexChrs_noFilters_KoryAnalysis_21Kgenes.txt") %>% 
      as_tibble(rownames = "ensembl_gene_id") %>% 
      rename_at(2:ncol(.), ~substr(.x, start = 2, stop = nchar(.x))) %>% 
      clean_names() %>% 
      dplyr::rename_at(2:ncol(.), ~str_remove(.x, "x"))

    # COVARIATES
    df_covariates_all <- read.table("data/covariates/185samples_multinomialDx_46covariates_10evecs.txt", 
               header = TRUE,
               fill = TRUE,
               sep = "\t") %>% 
      as_tibble() %>% 
      clean_names() %>% 
      rename("sample" = "x") %>% 
      mutate_if(is.character, ~tolower(.x))
      
      # remove outlier!!




# reformat covariates -----------------------------------------------------


    # CLEAN UP COVARIATE DATA - COMBINE COVARIATE CATAEGORIES
    df_covariates_all_clean <- df_covariates_all %>% 
      
      dplyr::select(-contains(c("evec", "sq", "_1", "_2", "hallucinogens"))) %>% # no one was doing hallucinogens
      select(-etoh_results) %>% 
      
      dplyr::rename("anti_depressant" = "antidepressants") %>% 
      select(-anti_depress) %>% 
      
      dplyr::rename("anti_anxiety" = "sedative_hypnotic_anxiolitics") %>% 
      select(-benzos) %>% 
      
      dplyr::rename("cannabis" = "cannabinoids") %>% 
      select(-thc) %>% 
      
      mutate(mood_stabilizers = case_when(
        
        mood_stabilizers == "negative" & mood_stab == "negative" ~ "negative",
        mood_stabilizers == "positive" & mood_stab == "positive" ~ "positive",
        mood_stabilizers == "negative" & mood_stab == "" ~ "negative",
        mood_stabilizers == "positive" & mood_stab == "" ~ "positive",
        mood_stabilizers == "" & mood_stab == "negative" ~ "negative",
        mood_stabilizers == "" & mood_stab == "positive" ~ "positive",
        
        TRUE ~ ""
        
      )) %>% select(-mood_stab) %>% 
      
      dplyr::rename("nicotine" = "nic_cot") %>% 
      select(-opiates)
    
    
    # CONVERT TO NUMERIC
    pos_neg_col_range <- 
      c(seq(from = which(colnames(df_covariates_all_clean) == "nicotine") - 1,
            to = which(colnames(df_covariates_all_clean) == "cocaine") - 1),
        seq(from = which(colnames(df_covariates_all_clean) == "alcohol") - 1,
            to = which(colnames(df_covariates_all_clean) == "non_psychiatric") - 1)) 
    
    df_covariates_all_clean_numeric <- df_covariates_all_clean %>% 
      select(-new_dx) %>% 
      mutate(gender = ifelse(gender == "f", 1, 0)) %>% 
      mutate(race = ifelse(race == "aa", 1, 0)) %>% 
      mutate(manner = recode(manner, 
                             accident = 1, 
                             homicide = 2, 
                             natural = 3, 
                             suicide = 4, 
                             undetermined = NA_real_)) %>% 
      mutate(smoker = ifelse(smoker == "yes", 1, 0)) %>% 
      mutate(source = recode(source,
                             `dc meo` = 1,
                             `northern virginia meo` = 2,
                             `central virginia meo` = 3,
                             `funeral home donation` = 3,
                             `other meo` = 3)) %>% 
      mutate(marital_status = recode(marital_status, 
                                     divorced = 1, 
                                     married = 2, 
                                     separated = 3, 
                                     single = 4, 
                                     widowed = 5,
                                     undetermined = NA_real_)) %>% 
      mutate_at(pos_neg_col_range, ~recode(.x, positive = 1, negative = 0)) %>% 
      mutate(suicide = ifelse(manner == 4, 1, 0))
    
    # SAVE
    # save(df_covariates_all_clean, df_covariates_all_clean_numeric,
    #      file = "data/covariates/185_all_covariates_clean.Rdata")
    
    
# raw count PCA, all genes ---------------------------------------------------------------------


    # RUN PCA
    pca <- prcomp(df_raw_counts %>% 
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
    
    # PLOT
    # p_pca <- 
      
    df_pca %>% 
      ggplot(aes(x = pc1, y = pc2)) +
      geom_point(color = "black", size = 2.75) +
      geom_point(aes(color = g_cpercent, text = sample), size = 2) +
      
      # label samples
      # geom_text_repel(aes(label = sample)) +
      
      # color gradient for continuous color variables
      scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
      
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
    
    # ggplotly(p_pca, tooltip = "text")
    
    
    # CORRELATE EACH PC WITH EACH COVARIATE
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
    

    # FIND SIGNIFICANT CORRELATIONS
    df_pc_covar_cor_sig <- df_pc_covar_cor %>% 
      mutate(pc = str_remove(pc, "pc") %>% as.numeric()) %>% 
      left_join(df_var %>% 
                  pivot_longer(2:ncol(.), names_to = "metric", values_to = "value") %>% 
                  filter(metric == "proportion_of_variance"), by = "pc") %>% 
      dplyr::select(-metric) %>% 
      relocate(value, .before = covariate) %>% 
      dplyr::rename("pc_perc_var" = "value") %>% 
      mutate(pc_perc_var = 100*pc_perc_var) %>% 
      mutate(p_adj = p.adjust(p_val, method = "BH")) %>% 
      filter(p_adj < 0.05 & pc_perc_var > 0.1)
    
    # df_pc_covar_cor_sig %>% 
    #   write.csv("outputs/20220609_pc_14k_genes_vsd_regress_RIN_race_GC_libraryBatch_significant_covariate_correlations_ALL_COV_NO_OUTLIER.csv")
      
    
    
    
# count normalized, vsd transformed PCA, all genes ---------------------------------------------------------------------
    
    
    # GET VSD COUNTS
    
        # convert raw counts to matrix
        m_raw_counts <- df_raw_counts %>% 
          
          # convert to numeric matrix
          as.data.frame %>% 
          column_to_rownames("ensembl_gene_id") %>% 
          as.matrix
        
        # check that rows and columns are aligned
        all(df_covariates_all_clean$sample == colnames(m_raw_counts))
        
        # create deseq object
        dds <- DESeqDataSetFromMatrix(countData = m_raw_counts,
                                      colData = df_covariates_all_clean,
                                      design = ~ new_dx)
        
        
        # estimate DESeq size (normalization) factors
        dds <- estimateSizeFactors(dds)
        
        # extraction
        m_dds_counts <- counts(dds, normalized = TRUE)
        
        # generate vsd
        vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
    
        # get vsd for each gene
        df_vsd <- vsd %>% assay %>% as.data.frame %>% rownames_to_column("ensembl_gene_id") %>% as_tibble()
        
    # RUN PCA
    pca <- prcomp(df_vsd %>% 
                    as.data.frame %>%
                    column_to_rownames("ensembl_gene_id"))
    
    pca <- prcomp(m_dds_counts)
    
    
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
    
    # PLOT
    # p_pca <- 
    
    df_pca %>% 
      ggplot(aes(x = pc1, y = pc2)) +
      geom_point(color = "black", size = 2.75) +
      geom_point(aes(color = rin_acsg, text = sample), size = 2) +
      
      # label samples
      # geom_text_repel(aes(label = sample)) +
      
      # color gradient for continuous color variables
      scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
      
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
    
    # ggplotly(p_pca, tooltip = "text")
    
    # 3D PLOT
    # pc1 <- df_pca$pc1
    # pc2 <- df_pca$pc2
    # pc3 <- df_pca$pc3
    # mapped_percent <- df_pca$mapped_percent
    # 
    # plot_ly(x = pc1, y = pc2, z = pc3, type = "scatter3d", mode = "markers", color = mapped_percent)
    
    
    # CORRELATE EACH PC WITH EACH COVARIATE
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
    
    
    # FIND SIGNIFICANT CORRELATIONS
    df_pc_covar_cor_sig <- df_pc_covar_cor %>% 
      mutate(pc = str_remove(pc, "pc") %>% as.numeric()) %>% 
      left_join(df_var %>% 
                  pivot_longer(2:ncol(.), names_to = "metric", values_to = "value") %>% 
                  filter(metric == "proportion_of_variance"), by = "pc") %>% 
      dplyr::select(-metric) %>% 
      relocate(value, .before = covariate) %>% 
      dplyr::rename("pc_perc_var" = "value") %>% 
      mutate(pc_perc_var = 100*pc_perc_var) %>% 
      mutate(p_adj = p.adjust(p_val, method = "BH")) %>% 
      filter(p_adj < 0.05 & pc_perc_var > 0.1)
    
    # df_pc_covar_cor_sig %>% 
    #   write.csv("outputs/20220609_pc_14k_genes_vsd_regress_RIN_race_GC_libraryBatch_significant_covariate_correlations_ALL_COV_NO_OUTLIER.csv")
    


    

# are age at death and drug use correlated with diagnosis? ----------------

    
    # look at missing data
    naniar::gg_miss_upset(df_covariates_all_clean_numeric)
    naniar::gg_miss_upset(df_covariates_all_clean_numeric %>% 
                    select(-c(anti_epileptics, anti_histamines, anticholinergics, nicotine, cannabis, non_psychiatric)))

    # CORRELATE COVARIATE MATRIX
    cov_cor <- df_covariates_all_clean_numeric %>% 
      as.data.frame %>% 
      column_to_rownames("sample") %>% 
      as.matrix %>% 
      Hmisc::rcorr(type = "spearman")

    df_cov_cor <- cov_cor$r %>% 
      as.data.frame %>% 
      rownames_to_column("covariate1") %>% 
      as_tibble() %>% 
      pivot_longer(2:ncol(.), names_to = "covariate2", values_to = "spearmans_r") %>% 
      
      left_join(
        
        cov_cor$P %>% 
          as.data.frame %>% 
          rownames_to_column("covariate1") %>% 
          as_tibble() %>% 
          pivot_longer(2:ncol(.), names_to = "covariate2", values_to = "p_value")
        
      ) %>% 
      
      filter(spearmans_r != 1) %>% 
      mutate(p_adj = p.adjust(p_value, method = "BH"))
    
    df_cov_cor %>% filter(p_adj < 0.05) %>% write.csv("outputs/20220609_184samples_covariate_correlations.csv")
    
    
    # PLOT
          technical_cov <- c("rn_aextraction_batch", "library_batch", "g_cpercent",
                             "five_prime_three_prime_bias", "rin_acsg", "age_death")
    
          df_cov_cor_zscore <- df_cov_cor %>% 
            mutate(z_score = (spearmans_r - mean(spearmans_r, na.rm = TRUE))/sd(spearmans_r, na.rm = TRUE)) %>% 
            filter(p_adj < 0.05) %>% 
            
            # remove technical covariates
            filter(!(covariate1 %in% technical_cov) & !(covariate2 %in% technical_cov))
          
          m_cov_cor_zscore <- df_cov_cor_zscore %>% 
            pivot_wider(id_cols = covariate1, names_from = covariate2, values_from = z_score) %>% 
            as.data.frame %>% 
            column_to_rownames("covariate1") %>% 
            as.matrix
          
          # cluster dendrogram
          cov_tree <- m_cov_cor_zscore %>% 
            dist %>% 
            hclust
          plot(cov_tree)
          cov_dendro <- as.dendrogram(cov_tree)
          labels <- labels(cov_dendro)
          labels <- labels[!labels %in% technical_cov]
          
          # plot heatmap
          df_cov_cor_zscore %>% 
            
            ggplot(aes(x = covariate1, y = covariate2)) +
            geom_tile(aes(fill = z_score)) +
            scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), na.value = "white") +
            scale_x_discrete(limits = labels) +
            scale_y_discrete(limits = labels) + 
            labs(x = "", y = "") +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.9),
              axis.text = element_text(size = 10))
    
          # plot multidx correlations
          dx_sig_cor <- df_cov_cor_zscore %>% filter(covariate1 == "multi_nomial_dx") %>% pull(covariate2)

          df_covariates_all_clean %>% 
            mutate_if(is.character, ~ifelse(.x == "", "NA", .x)) %>% 
            mutate(marital_status = ifelse(marital_status == "unknown", "NA", marital_status)) %>% 
            mutate(new_dx = factor(new_dx, levels = c("control", "bipolar", "mdd", "schizo"))) %>% 
            dplyr::select(new_dx, contains(dx_sig_cor)) %>% 
            pivot_longer(2:ncol(.), values_to = "value", names_to = "covariate") %>% 
            filter(value != "NA") %>% 
            
            # find percent for each covariate
            group_by(new_dx, covariate) %>% 
            count(value) %>% 
            mutate(percent = 100*n/sum(n)) %>% 
            
            # plot
            ggplot(aes(x = value, fill = new_dx)) +
            geom_bar(aes(y = percent), stat = "identity", position = "dodge") +
            facet_wrap(~covariate, scales = "free") +
            labs(x = "")
          
          glm(smoker ~ multi_nomial_dx, 
              family = "binomial", 
              df_covariates_all_clean_numeric %>% 
                mutate(multi_nomial_dx = factor(multi_nomial_dx, levels = c("1", "2", "3", "4")))) %>% 
            tidy(exp = TRUE)
          
          # plot race correlations
          race_sig_cor <- df_cov_cor_zscore %>% filter(covariate1 == "race") %>% pull(covariate2)

          df_covariates_all_clean_numeric %>% 
            mutate_if(is.character, ~ifelse(.x == "", "NA", .x)) %>% 
            mutate(race = factor(race, levels = c("1", "0"))) %>% 
            dplyr::select(race, contains(race_sig_cor)) %>% 
            pivot_longer(2:ncol(.), values_to = "value", names_to = "covariate") %>% 
            filter(value != "NA") %>% 
            
            # find proportion count for each covariate
            group_by(race, covariate) %>% 
            count(value) %>% 
            mutate(percent = 100*n/sum(n)) %>% 
            
            ggplot(aes(x = value, fill = race)) +
            geom_bar(aes(y = percent), stat = "identity", position = "dodge") +
            facet_wrap(~covariate, scales = "free") +
            labs(x = "", y = "percent")
          
          
          
    
# correlation of gc_percent and mapped_percent ----------------------------
          
          
          mp_gc_r <- paste0("r = ", 
                            df_cov_cor %>% filter(covariate1 == "mapped_percent", covariate2 == "g_cpercent") %>% pull(spearmans_r) %>% round(3))
          mp_gc_p <- paste0("p = ", 
                            df_cov_cor %>% filter(covariate1 == "mapped_percent", covariate2 == "g_cpercent") %>% pull(p_adj) %>% round(3))
          
          mp_rin_r <- paste0("r = ", 
                            df_cov_cor %>% filter(covariate1 == "mapped_percent", covariate2 == "rin_acsg") %>% pull(spearmans_r) %>% round(3))
          mp_rin_p <- paste0("p = ", 
                            df_cov_cor %>% filter(covariate1 == "mapped_percent", covariate2 == "rin_acsg") %>% pull(p_adj) %>% round(3))
          
          rin_gc_r <- paste0("r = ", 
                            df_cov_cor %>% filter(covariate1 == "rin_acsg", covariate2 == "g_cpercent") %>% pull(spearmans_r) %>% round(3))
          rin_gc_p <- paste0("p = ", 
                            df_cov_cor %>% filter(covariate1 == "rin_acsg", covariate2 == "g_cpercent") %>% pull(p_adj) %>% round(3))
          
          df_covariates_all_clean_numeric %>% 
            ggplot(aes(x = rin_acsg, y  = g_cpercent)) +
            geom_point(size = 3, shape = 1) +
            geom_smooth(method = "lm") +
            annotate(geom = "text", x = 8.5, y = 52, label = paste0(rin_gc_r, "; ", rin_gc_p), size = 8)
          
          
    

      
    
    
    