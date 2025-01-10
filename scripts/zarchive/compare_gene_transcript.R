


# libraries ---------------------------------------------------------------


  library(tidyverse)
  library(janitor)
  library(ggtranscript)


# set theme for plots -----------------------------------------------------

theme_set(theme_bw() +
            theme(plot.title = element_text(size = 18),
                  axis.title = element_text(size = 15),
                  axis.text = element_text(size = 12),
                  strip.text = element_text(size = 15),
                  legend.title = element_text(size = 15),
                  legend.text = element_text(size = 12)))


# data --------------------------------------------------------------------


  base_dir <- "/Users/work/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/projects/WGCNA/acsg_wgcna/"
  prefix_gene <- "20230228_185samples_19kgenes_vst_qSVA123567_MP_RNAbatch_Race_resids"
  prefix_transcript <- "20230315_185samples_25kRAREtranscripts_vst_qSVA123567_MP_RNAbatch_Race_resids"
  soft_power <- 3

# LOAD GENE OBJECTS  
  load(paste0("objects/", prefix_gene, ".RDS")) # df_vsd_regress
  load(paste0("objects/", prefix_gene, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

  df_modules_gene <- df_modules %>% 
    filter(min_size == 40 & cut_height == 0.97) %>% 
    unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module)
  
# LOAD RARE TRANSCRIPTS
  load(paste0(base_dir, "objects/", prefix_transcript, ".RDS")) # df_vsd_regress_rare
  load(paste0(base_dir, "objects/", prefix_transcript, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules
  
  df_modules_transcript <- df_modules %>% 
    filter(min_size == 30 & cut_height == 0.9986) %>% 
    unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module)
  
# LOAD OTHER INFO
  load("objects/transcript_gene_go_term.RDS") # df_transcript_gene_go_term
  load("objects/ensembl_id_to_gene_symbol.Rdata") # df_ensembl_to_symbol

# LOAD HUMAN GENOME
  load("objects/df_hsapiens_genome.RDS") # df_hsapiens_genome
  
  
# load CCA res ------------------------------------------------------------

  

# TRANSCRIPT
  cca_dir_transcript <- paste0(base_dir, "RCCA_toolkit/s185_25kRAREtranscripts_covarOpioidsAnxChol_regressAgedeath/")
  model_mat_transcript <- read_mat(paste0(cca_dir_transcript, "framework/cca_holdout1-0.20_subsamp5-0.20/res/level2/model_1.mat"))
  
  transcript_universe <- df_vsd_regress_rare %>%     
    dplyr::select(ensembl_transcript_id, sample, resids) %>% 
    left_join(df_modules_transcript) %>% 
    arrange(module) %>% 
    pivot_wider(id_cols = sample, names_from = ensembl_transcript_id, values_from = resids) %>%
    dplyr::select(-sample) %>% 
    colnames()
  
  df_x_res_transcript <- tibble(ensembl_transcript_id = transcript_universe,
                                weight = model_mat_transcript$wX) %>% 
    left_join(df_hsapiens_genome %>% 
                dplyr::select(transcript_id, transcript_name,
                              gene_id, gene_name) %>% 
                dplyr::rename("ensembl_transcript_id" = "transcript_id") %>% 
                distinct()) %>% 
    left_join(df_modules_transcript) %>% 
    mutate(abs_weight = abs(weight)) %>% 
    arrange(-abs_weight) %>% 
    mutate(ensembl_transcript_id = factor(ensembl_transcript_id, 
                                          levels = unique(.$ensembl_transcript_id))) 
  
# GENE
  cca_dir_gene <- paste0(base_dir, "RCCA_toolkit/s185_19kgenes_covarOpioidsAnx_regressBrainweightAgedeath/")
  model_mat_gene <- read_mat(paste0(cca_dir_gene, "framework/cca_holdout1-0.20_subsamp5-0.20/res/level1/model_1.mat"))
 
  gene_universe <- df_vsd_regress %>%     
    dplyr::select(ensembl_gene_id, sample, resids) %>% 
    left_join(df_modules_gene) %>% 
    arrange(module) %>% 
    pivot_wider(id_cols = sample, names_from = ensembl_gene_id, values_from = resids) %>%
    dplyr::select(-sample) %>% 
    colnames()
    
  df_x_res_gene <- tibble(ensembl_gene_id = gene_universe,
                          weight = model_mat_gene$wX) %>% 
    left_join(df_modules_gene) %>% 
    mutate(abs_weight = abs(weight)) %>% 
    arrange(-abs_weight) %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = unique(.$ensembl_gene_id)))
  
  
# What modules do the rare transcripts end up in? -------------------------
  
  
# PERCENT OF MODULES THAT ARE RARE TRANSCRIPTS
  df_compare_gene_tx <- df_modules_transcript %>% 
    left_join(df_transcript_gene_go_term) %>% 
    dplyr::rename("gene_symbol" = "external_gene_name",
                  "transcript_module" = "module",
                  "transcript_color" = "color") %>% 
    dplyr::select(-go_id, -term, -mod_set) %>% 
    left_join(df_modules_gene %>% dplyr::select(-mod_set))  %>% 
    dplyr::rename("gene_module" = "module",
                  "gene_color" = "color") %>% 
    distinct()
  
  df_compare_gene_tx %>% 
    filter(gene_module == 0) %>% 
    count(transcript_module)
  
  

# compare to Akula et al., 2021 results -----------------------------------

# LOAD NIRMALA'S DE RESULTS  
  load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res  
  

# TOP 1% WEIGHTED CCA TRANSCRIPTS
  top_cca_transcripts <- df_x_res_transcript %>% 
    head(0.01*nrow(.)) %>% 
    pull(ensembl_transcript_id) %>% 
    as.character()
  
# TOP 1% WEIGHT CCA GENES 
  top_cca_genes <- df_x_res_gene %>% 
    head(0.01*nrow(.)) %>% 
    pull(ensembl_gene_id) %>% 
    as.character()
  
df_akula_res %>% 
  filter(id %in% c(top_cca_transcripts, top_cca_genes) & padj < 0.05) %>% 
  left_join(df_x_res_gene %>% 
              dplyr::rename("id" = "ensembl_gene_id") %>% 
              dplyr::select(id, weight, module) %>% 
              mutate(module = as.character(module))
            ) %>% 
  left_join(df_x_res_transcript %>% 
              mutate(module = as.character(module)) %>% 
              dplyr::rename("id" = "ensembl_transcript_id",
                            "cca_weight" = "weight",
                            "wgcna_module" = "module") %>% 
              dplyr::select(id, cca_weight, wgcna_module)
            ) %>% 
  mutate(cca_weight = ifelse(is.na(cca_weight), weight, cca_weight),
         wgcna_module = ifelse(is.na(wgcna_module), module, wgcna_module)) %>% 
  dplyr::select(-weight, -module) #%>% 

  
  write.csv(paste0(base_dir, "outputs/cca_overlap_with_akula_DE_res.csv"),
            row.names = FALSE)
  
  
  
# CCA ---------------------------------------------------------------------


  top_cca_transcript_to_gene <- df_hsapiens_genome %>% 
    filter(transcript_id %in% top_cca_transcripts) %>% 
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>% 
    pull(gene_name) %>% 
    unique
  
  df_x_res_gene %>% 
    head(0.01*nrow(.)) %>% 
    
    filter(gene_symbol %in% top_cca_transcript_to_gene)

# PLOT
  # specify transcripts
  transcripts <- df_x_res_transcript %>%
    filter(gene_name == "EPS15") %>%
    pull(ensembl_transcript_id) %>%
    as.character()
  
  # extract exons
  exons <- df_hsapiens_genome %>% 
    dplyr::filter(transcript_id %in% transcripts, type == "exon") %>% 
    mutate(transcript_name = ifelse(is.na(transcript_name),
                                    transcript_id,
                                    transcript_name))
  
  top_transcript_name <- exons %>% 
    filter(transcript_id ==
             (df_x_res_transcript %>%
                head(0.01*nrow(.)) %>%
                filter(gene_name == "EPS15") %>%
                pull(ensembl_transcript_id) %>%
                as.character())
    ) %>%
    pull(transcript_name) %>% 
    unique
  
  exons %>% 
    mutate(transcript_name = ifelse(transcript_name %in% top_transcript_name,
                                    paste0("*", transcript_name),
                                    transcript_name)) %>% 
    
    ggplot(aes(xstart = start, xend = end, y = transcript_name)) +
    geom_range(aes(fill = transcript_biotype)) +
    geom_intron(data = to_intron(exons %>% 
                                   mutate(transcript_name = 
                                            ifelse(transcript_name %in% top_transcript_name,
                                                   paste0("*", transcript_name),
                                                   transcript_name)),
                                 "transcript_name"),
                aes(strand = strand)) +
    facet_wrap(vars(gene_name), scales = "free") +
    ylab("") +
    theme(legend.position = "bottom")
  
# How many of the top gene hits have rare transcripts?  
  top_cca_gene <- df_x_res_gene %>% 
    head(0.01*nrow(.)) %>% 
    pull(ensembl_gene_id) %>% 
    as.character()
  
  top_cca_gene_to_transcript <- df_hsapiens_genome %>% 
    filter(gene_id %in% top_cca_gene) %>% 
    pull(transcript_id) %>% 
    unique()
  
  df_x_res_transcript %>% 
    head(0.1*nrow(.)) %>% 
    filter(ensembl_transcript_id %in% top_cca_gene_to_transcript)
    
  
# PLOT
  
  TOIs <- c( "ENST00000498255",
             "ENST00000594077",
             "ENST00000467906",
             "ENST00000490310"
  )
  
  corresponding_genes <- df_hsapiens_genome %>% 
    filter(transcript_id %in% TOIs) %>% 
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>% 
    pull(gene_name) %>% unique
  
  # specify transcripts
  
  # transcripts <- df_x_res_transcript %>% 
  #   filter(gene_symbol == "DGKH") %>% 
  #   pull(ensembl_transcript_id) %>% 
  #   as.character()
  
  transcripts <- df_x_res_transcript %>% 
    filter(gene_name %in% corresponding_genes |
             gene_id %in% corresponding_genes) %>% 
    pull(ensembl_transcript_id) %>% 
    as.character()
  
  # extract exons
  exons <- df_hsapiens_genome %>% 
    dplyr::filter(transcript_id %in% transcripts, type == "exon") %>% 
    mutate(transcript_name = ifelse(is.na(transcript_name),
                                    transcript_id,
                                    transcript_name))
  
  top_transcript_name <- exons %>% 
    filter(transcript_id %in% TOIs) %>% 
    # filter(transcript_id == 
    #          (df_x_res_transcript %>% 
    #             head(0.01*nrow(.)) %>%
    #             filter(gene_symbol == "DGKH") %>%
    #             pull(ensembl_transcript_id) %>%
    #             as.character())
    # ) %>% 
    pull(transcript_name) %>% 
    unique
  
  exons %>% 
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>% 
    mutate(transcript_name = ifelse(transcript_name %in% top_transcript_name,
                                    paste0("*", transcript_name),
                                    transcript_name)) %>% 
    
    ggplot(aes(xstart = start, xend = end, y = transcript_name)) +
    geom_range(aes(fill = transcript_biotype)) +
    geom_intron(data = to_intron(exons %>% 
                                   mutate(transcript_name = 
                                            ifelse(transcript_name %in% top_transcript_name,
                                                   paste0("*", transcript_name),
                                                   transcript_name)),
                                 "transcript_name"),
                aes(strand = strand)) +
    facet_wrap(vars(gene_name), scales = "free") +
    ylab("") +
    theme(legend.position = "bottom")
  
  
  

  
