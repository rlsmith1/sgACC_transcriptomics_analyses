
########################################################################################

# GENE LEVEL: run GO and cell-type enrichment on WGCNA modules

########################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/supplement/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WGCNA_res.R"))


# Run GO on each WGCNA module --------------------------------------------

  
## Run GO pathway overrepresentation analysis on each module
modules <- df_modules_filt %>% pull(module) %>% unique
doParallel::registerDoParallel()
df_mods_go <- map_dfr(
  .x = modules,
  .f = ~ enrichGO(gene = df_modules_filt %>% dplyr::filter(module == .x) %>% pull(ensembl_gene_id),
                  OrgDb = "org.Hs.eg.db",
                  universe = ensembl_gene_ids,
                  keyType = "ENSEMBL",
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  readable = TRUE
  ) %>% 
    as_tibble() %>% 
    mutate(module = .x, .before = 1)
)


## Save results for figures & tables ##


# Figures
save(df_mods_go, file = paste0(analysis_objects_dir, "WGCNA_module_GO.RDS"))

# Downstream analyses
save(df_mods_go, file = paste0(base_dir, "objects/", prefix, sep = "_", "SIGNED_SFT", soft_power, "_GO_RES.RDS"))
  



# Run topic modeling on GO results for each module ------------------------------------------------------

#load(paste0(analysis_objects_dir, "WGCNA_module_GO.RDS"))

## Unnest tokens to get GO terms within each module (group)
df_go_tokens <- df_mods_go %>% 
    dplyr::filter(pvalue < 0.05) %>% 
    dplyr::filter(!is.na(Description)) %>% 
    group_by(module) %>% 
    unnest_tokens(word, Description)

## Create a sparse matrix using only meaningful words that appear at least twice within module
vague_terms <- c("positive", "negative", "regulation", "response", "protein", "activity",
                 "pathway", "process", "involved", "signaling")
go_sparse <- df_go_tokens %>%
    dplyr::count(module, word) %>%
    mutate(module = as.numeric(as.character(module %>% str_remove("geneM")))) %>% 
    anti_join(stop_words) %>% 
    dplyr::filter(n > 2 & !(word %in% vague_terms)) %>%
    cast_sparse(module, word, n)

dim(go_sparse)

## Fit topic model
set.seed(456)
n_topics <- 5
topic_model <- stm(go_sparse, K = n_topics, verbose = FALSE)


## Save the topic model and sparse matrix for plotting
save(go_sparse, topic_model,
     file = paste0(analysis_objects_dir, "topic_modeling_res.Rdata"))



# Cell-type enrichment ----------------------------------------------------

## Run hypergeometric test of each module with each Lake cell type
df_mods_cell_type_hypergeometric <- expand_grid(
    cell_type = names(cell_types),
    module = modules
) %>% 
    mutate(
        hypergeometric_res = map2(
            .x = cell_type,
            .y = module,
            .f = ~ f_module_hypergeometric(
                gene_list = cell_types[[.x]], 
                gene_list_name = .x, 
                wgcna_module = .y, 
                gene_universe = ensembl_gene_ids, 
                module_data = df_modules_filt
            ) %>% 
                dplyr::select(-module)
        )
    ) %>% 
    unnest(cols = c(hypergeometric_res)) %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr"))

## Save for plotting
save(df_mods_cell_type_hypergeometric,
     file = paste0(analysis_objects_dir, "cell_type_res.RDS"))



# Approximate developmental trajectory of each module using PsychENCODE data ----------------------------------------------

load(paste0(base_dir, "objects/psychencode_development_expr_data.Rdata")) # df_psychencode_metadata, df_psychencode_expr

## Identify frontal cortex samples
df_sample_window <-  df_psychencode_metadata %>% 
    dplyr::filter(region_broad == "frontal_cortex")
frontal_cortex_samples <- df_sample_window %>% 
    pull(id)

## Filter expression data for samples and genes of interest  
df_psychencode_expr_filt <- df_psychencode_expr %>% 
    dplyr::filter(GENE %in% ensembl_gene_ids) %>% 
    dplyr::select(GENE, all_of(frontal_cortex_samples)) %>% 
    dplyr::rename("ensembl_gene_id" = "GENE")

## Combine expresion data and module data
df_expr_mod <- df_psychencode_expr_filt %>% 
    pivot_longer(starts_with("HS"), names_to = "id", values_to = "expression") %>% 
    left_join(df_sample_window) %>% 
    left_join(df_modules_filt) %>% 
    dplyr::select(ensembl_gene_id, module, id, window, expression) 

## Find median expression value of each sample for each module
df_expr_mod_median <- df_expr_mod %>% 
    group_by(module, id, window) %>%
    summarise(median_expr = median(expression))

## Save for plotting
save(df_expr_mod_median,
     file = paste0(analysis_objects_dir, "developmental_trajectory.RDS"))

