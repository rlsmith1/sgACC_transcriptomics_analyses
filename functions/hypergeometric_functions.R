
# General hypergeometric function -----------------------------------------

#=====================================================================================================================#

### GOAL: Test hypergeometric overlap between two gene lists

### INPUTS ###

## gene_list1 = list of genes (can be any type of ID, e.g., Ensembl ID or HGNC symbol)
## gene_list2 = list of genes (must be same type as gene_list1)
## gene_universe = list of all the genes that were included in the analysis (same ID type as gene_list1 & gene_list2)

### OUTPUTS ###

## overlap_genes = list of genes that were in both gene_list1 and gene_list2
## p_value = hypergeometric p-value (significant of overlap given total pool of genes)

#=====================================================================================================================#

f_hypergeometric <- function(gene_list1, gene_list2, gene_universe) {
    
    # Include only genes that were part of universe
    gene_list1_in_universe <- intersect(gene_universe, gene_list1) %>% unique
    gene_list2_in_universe <- intersect(gene_universe, gene_list2) %>% unique
    
    # Identify the overlap between the published gene list (including only those in our universe) and the gene list of interest
    overlap_genes <- intersect(gene_list1_in_universe, gene_list2_in_universe)
    
    # Run hypergeometric test
    phyper <- 1 - phyper(q = length(overlap_genes), 
                         m = length(gene_list2_in_universe), 
                         n = length(gene_universe) - length(gene_list2_in_universe), 
                         k = length(gene_list1_in_universe)
    ) 
    
    # Return a list of common genes and the hypergeometric p-value
    return(list("overlap_genes" = overlap_genes, "p_value" = phyper))
    
}



# WGCNA module hypergeometric ----------------------------------------

#=====================================================================================================================#

### GOAL: Specific case of hypergeometric function to test overlap of a gene list with WGCNA module

### INPUTS ###

## gene_list = list of genes (can be any type of ID, e.g., Ensembl ID or HGNC symbol)
## gene_list_name = name of the gene list that you want to include for the output table
## wgcna_module = name of the WGCNA-derived module of interest
## gene_universe = list of all the genes that were included in the analysis (same ID type as gene_list1 & gene_list2)
## module_data = a dataframe or tibble of gene module assignments. Must have a column called "module" (with WGCNA module names) and a column called "ensembl_gene_id" (with gene IDs)

### OUTPUTS ###

## a tibble with 5 columns
    # gene_list = name of the gene_list (identified by gene_list_name)
    # module = name of the WGCNA module (wgcna_module)
    # p_value = hypergeometric p-value
    # total_n = total number of genes in the WGCNA module
    # overlap_n = number of genes in the overlap (both in gene_list and wgcna_module)

#=====================================================================================================================#

## Write function to test the hypergeometric overlap of a given gene list with each WGCNA module
f_module_hypergeometric <- function(gene_list, gene_list_name, wgcna_module, gene_universe, module_data = df_modules_filt) {
    
    # Filter module tibble for module of interest
    module_genes <- module_data %>% dplyr::filter(module == wgcna_module) %>% pull(ensembl_gene_id)
    
    # Run hypergeometric test on gene list & module
    hypergeometric_res <- f_hypergeometric(gene_list1 = module_genes, 
                                           gene_list2 = gene_list,
                                           gene_universe = gene_universe
    )
    
    # Create tibble of results to output
    df_res <- tibble(
        gene_list = gene_list_name,
        module = wgcna_module,
        p_value = hypergeometric_res$p_value,
        total_n = length(module_genes),
        overlap_n = length(hypergeometric_res$overlap_genes)
    )
    return(df_res)
    
}


# WGCNA module hypergeometric (TRANSCRIPT CASE) ----------------------------------------

#=====================================================================================================================#

### GOAL: Specific case of hypergeometric function to test overlap of a gene list with WTCNA module (transcript case)

## Same inputs and outputs as above but for specific case of transcript use (rather than genes)

#=====================================================================================================================#

## Write function to test the hypergeometric overlap of a given gene list with each WGCNA module
f_module_hypergeometric_transcript <- function(transcript_list, transcript_list_name, wgcna_module, transcript_universe = ensembl_ids) {
    
    # Filter module tibble for module of interest
    module_transcripts <- df_modules_filt %>% dplyr::filter(module == wgcna_module) %>% pull(ensembl_transcript_id)
    
    # Run hypergeometric test on transcript list & module
    hypergeometric_res <- f_hypergeometric(gene_list1 = module_transcripts, 
                                           gene_list2 = transcript_list,
                                           gene_universe = transcript_universe
    )
    
    # Create tibble of results to output
    df_res <- tibble(
        gene_list = transcript_list_name,
        module = wgcna_module,
        p_value = hypergeometric_res$p_value,
        total_n = length(module_transcripts),
        overlap_n = length(hypergeometric_res$overlap_genes)
    )
    return(df_res)
    
}

