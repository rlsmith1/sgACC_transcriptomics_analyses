
library(tidyverse)
library(readxl)
library(janitor)

df_akula_res <- 
  
  # gene level
  readxl::read_excel("data/Akula_etal_2021_DEGs.xlsx", sheet = 4) %>% 
  dplyr::select(1:6) %>% 
  row_to_names(1) %>%
  mutate(level = "gene",
         comparison = "schizo_ctrl",
         .before = 1) %>% 
  clean_names %>% 
  dplyr::rename("id" = "ensemble_gene_id",
                "symbol" = "hgnc_symbol") %>% 
  
  bind_rows(
    
    readxl::read_excel("data/Akula_etal_2021_DEGs.xlsx", sheet = 5) %>% 
      dplyr::select(1:6) %>% 
      row_to_names(1) %>%
      mutate(level = "gene",
             comparison = "bipolar_ctrl",
             .before = 1) %>% 
      clean_names %>% 
      dplyr::rename("id" = "ensemble_gene_id",
                    "symbol" = "hgnc_symbol")
    
  ) %>% 
  
  bind_rows(
    
    readxl::read_excel("data/Akula_etal_2021_DEGs.xlsx", sheet = 6) %>% 
      dplyr::select(1:6) %>% 
      row_to_names(1) %>%
      mutate(level = "gene",
             comparison = "mdd_ctrl",
             .before = 1) %>% 
      clean_names %>% 
      dplyr::rename("id" = "ensemble_gene_id",
                    "symbol" = "hgnc_symbol")
    
  ) %>% 
  
  bind_rows(
    
    readxl::read_excel("data/Akula_etal_2021_DEGs.xlsx", sheet = 7) %>% 
      dplyr::select(1:6) %>% 
      row_to_names(1) %>%
      mutate(level = "gene",
             comparison = "schizo_bipolar",
             .before = 1) %>% 
      clean_names %>% 
      dplyr::rename("id" = "ensemble_gene_id",
                    "symbol" = "hgnc_symbol")
    
  )  %>% 
  
  bind_rows(
    
    readxl::read_excel("data/Akula_etal_2021_DEGs.xlsx", sheet = 8) %>% 
      dplyr::select(1:6) %>% 
      row_to_names(1) %>%
      mutate(level = "gene",
             comparison = "schizo_mdd",
             .before = 1) %>% 
      clean_names %>% 
      dplyr::rename("id" = "ensemble_gene_id",
                    "symbol" = "hgnc_symbol")
    
  )  %>% 
  
  bind_rows(
    
    readxl::read_excel("data/Akula_etal_2021_DEGs.xlsx", sheet = 9) %>% 
      dplyr::select(1:6) %>% 
      row_to_names(1) %>%
      mutate(level = "gene",
             comparison = "bipolar_mdd",
             .before = 1) %>% 
      clean_names %>% 
      dplyr::rename("id" = "ensemble_gene_id",
                    "symbol" = "hgnc_symbol")
    
  )  %>% 
  
  # transcript level
  bind_rows(
    
    readxl::read_excel("data/Akula_etal_2021_DEGs.xlsx", sheet = 15) %>% 
      dplyr::select(1:6) %>% 
      row_to_names(1) %>%
      mutate(level = "transcript",
             comparison = "schizo_ctrl",
             .before = 1) %>% 
      clean_names %>% 
      dplyr::rename("id" = "transcript_name",
                    "symbol" = "hgnc_transcript_name_id")
    
  ) %>% 
  
  bind_rows(
    
    readxl::read_excel("data/Akula_etal_2021_DEGs.xlsx", sheet = 16) %>% 
      dplyr::select(1:6) %>% 
      row_to_names(1) %>%
      mutate(level = "transcript",
             comparison = "bipolar_ctrl",
             .before = 1) %>% 
      clean_names %>% 
      dplyr::rename("id" = "transcript_name",
                    "symbol" = "hgnc_transcript_name_id")
    
  )  %>% 
  
  bind_rows(
    
    readxl::read_excel("data/Akula_etal_2021_DEGs.xlsx", sheet = 17) %>% 
      dplyr::select(1:6) %>% 
      row_to_names(1) %>%
      mutate(level = "transcript",
             comparison = "mdd_ctrl",
             .before = 1) %>% 
      clean_names %>% 
      dplyr::rename("id" = "transcript_name",
                    "symbol" = "hgnc_transcript_name_id")
    
  ) %>% 
  
  # format
  mutate_at(vars(base_mean, log2fold_change, pvalue, padj), ~ as.numeric(.x))

save(df_akula_res, file = "objects/akula_et_al_DE_results.RDS") # df_akula_res

