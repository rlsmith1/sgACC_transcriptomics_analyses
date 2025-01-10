f_cell_type_overlap <- function(df,
                                module_no, 
                                cell_type) {
  
  # pull genes in module
  module_genes <- df %>% 
    dplyr::filter(module == module_no) %>% 
    pull(ensembl_gene_id)
  
  # pull genes in cell_type
  cell_type_genes <- df %>% 
    dplyr::filter(type == cell_type) %>% 
    pull(ensembl_gene_id)
  
  # find intersect of module with gene_universe_intersect genes
  module_intersect <- intersect(gene_universe_intersect, module_genes)
  
  # find intersect of cell_type with gene_universe_intersect genes
  cell_type_intersect <- intersect(gene_universe_intersect, cell_type_genes)
  
  # find intersect of module intersects
  module_cell_type_intersect <- intersect(module_intersect, cell_type_intersect)
  
  # set variables for hypergeometric function
  q <- (module_cell_type_intersect %>% length) - 1 # overlap between two modules - 1
  m <- module_intersect %>% length # module 1 
  n <- (gene_universe_intersect %>% length) - m # pop size - module 1
  k <- cell_type_intersect %>% length # module 2
  
  # return tibble
  tibble(module = module_no,
         type = cell_type,
         q = q,
         m = m,
         n = n,
         k = k,
         p_val = phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE))
  
}
