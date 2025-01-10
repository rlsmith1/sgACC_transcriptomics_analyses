


# libraries ---------------------------------------------------------------

  library(tidyverse)
  library(WGCNA)


# data --------------------------------------------------------------------


  prefix <- "20221209_185samples_19kgenes_vsd_qSVA17_resids"
  
  load(paste0("objects/", prefix, ".RDS")) # df_vsd_regress

  
# CONVERT TO MATRIX  
  m_vsd_regress <- df_vsd_regress %>% 
    pivot_wider(id_cols = sample,
                names_from = ensembl_gene_id,
                values_from = vsd) %>% 
    as.data.frame %>% 
    column_to_rownames("sample") 

  

# k-means elbow method ----------------------------------------------------

# Elbow method: percentage of variance explained as a function of the number of clusters
  
  set.seed(123)
  k_seq <- seq(2, 40, 1)
  
  df_kmeans <- map(
    
    .x = k_seq,
    .f = ~ tibble(
      
      k = .x,
      tot_within_ss = kmeans(m_vsd_regress, centers = .x)$tot.withinss
      
    )
    
  ) %>% bind_rows()
  
  p_elbow_vst <- df_kmeans %>% 
    ggplot(aes(x = k, y = tot_within_ss)) +
    geom_point() +
    ggtitle("Elbow method VST counts")

  # scaled data
  scaled_data <- scale(m_vsd_regress)

  df_kmeans_scaled <- map(
    
    .x = k_seq,
    .f = ~ tibble(
      
      k = .x,
      tot_within_ss = kmeans(scaled_data, centers = .x)$tot.withinss
      
    )
    
  ) %>% bind_rows()
  
  p_elbow_scaled <- df_kmeans_scaled %>% 
    ggplot(aes(x = k, y = tot_within_ss)) +
    geom_point() +
    ggtitle("Elbow method scaled data")
  
  # TOM
  adjacency <- adjacency(m_vsd_regress, power = 3)
  TOM <- TOMsimilarity(adjacency)
  diss_TOM <- 1 - TOM
  dist <- as.dist(diss_TOM)
  
  library(cluster)
  pam <- pam(diss_TOM, k = 2, diss = TRUE, variant = "faster")$silinfo$avg.width
  
  df_kmediods_tom <- map_dfr(
    
    .x = k_seq[1:29],
    .f = ~ tibble(
      
      k = .x,
      avg_width = pam(diss_TOM, 
                      k = .x, 
                      diss = TRUE, 
                      variant = "faster")$silinfo$avg.width
      
    )
    
  )
  
  p_elbow_tom <- df_kmediods_tom %>% 
    ggplot(aes(x = k, y = avg_width)) +
    geom_point() +
    ggtitle("Elbow method TOM")

p_elbow_vst + p_elbow_scaled + p_elbow_tom  
  
  

# k-means and Bayesian Information Criterion ------------------------------

adjacency <- adjacency(m_vsd_regress, power = 3)
TOM <- TOMsimilarity(adjacency)
diss_TOM <- 1 - TOM

save(diss_TOM, file = paste0("objects/", prefix, "_SFT3dissTOM.RDS"))

library(mclust)

doParallel::registerDoParallel()
clust <- Mclust(diss_TOM, G = 8, 
       modelNames = mclust.options("emModelNames"))
  
BIC <- clust$BIC
plot(BIC)
summary(BIC)

mod1 <- Mclust(diss_TOM, x = BIC)
summary(mod1, parameters = TRUE)
  
mod1 %>% str



# tidyclust ---------------------------------------------------------------



library(tidymodels)
library(tidyclust)

# data <- as.data.frame(diss_TOM)

data <- df_vsd_regress %>% 
  pivot_wider(id_cols = ensembl_gene_id,
              names_from = sample,
              values_from = vsd) %>% 
  as.data.frame %>% 
  column_to_rownames("ensembl_gene_id")

rec_spec <- recipe( ~ ., data = data) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_pca(all_numeric_predictors())

kmeans_spec <- k_means(num_clusters = tune())

wflow <- workflow() %>%
  add_recipe(rec_spec) %>%
  add_model(kmeans_spec)

grid <- tibble(num_clusters = 5:30)

set.seed(4400)
folds <- vfold_cv(data, v = 5)

gc()
res <- tune_cluster(
  wflow,
  resamples = folds,
  grid = grid,
  metrics = cluster_metric_set(sse_within_total, sse_total, sse_ratio)
)

res_metrics <- res %>% collect_metrics

res_metrics %>%
  filter(.metric == "sse_ratio") %>%
  ggplot(aes(x = num_clusters, y = mean)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  ylab("mean WSS/TSS ratio, over 5 folds") +
  xlab("Number of clusters") +
  scale_x_continuous(breaks = 1:30)


kclust <- kmeans(data, centers = 13)
    
df_kmeans_module <- kclust$cluster %>% 
  enframe %>% 
  dplyr::rename("ensembl_gene_id" = "name",
                "module" = "value") %>% 
  mutate(module = factor(module)) %>% 
  arrange(module)

df_kmeans_module %>% 
  count(module)


# SOURCE FUNCTION
source("functions/f_top_GO_modules.R")
load("objects/df_go_full_definitions.Rdata") # full go definitions

# PULL ALL WGCNA GENES AFTER FILTERING FOR BACKGROUND
l_gene_universe <- df_vsd_regress$ensembl_gene_id %>% unique

# CREATE GO ANNOTATION

# get gene info from ensembl database
hsapiens_ensembl <- useEnsembl(biomart = "genes",
                               dataset = "hsapiens_gene_ensembl")

# map each gene to GO term
m_go_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id"),
                  filters = "ensembl_gene_id",
                  values = l_gene_universe,
                  mart = hsapiens_ensembl)

# build annotation list
l_gene_2_GO <- m_go_ids[,c(1,3)] %>% unstack

library(janitor)
# RUN ENRICHMENT FUNCTION ON ALL MODULES
doParallel::registerDoParallel()

mods <- levels(df_kmeans_module$module)
df_mods_go <- map_dfr(
  
  .x =  mods,
  .f = ~ f_top_GO_modules(l_gene_module = df_kmeans_module %>%
                            filter(module == .x) %>%
                            pull(ensembl_gene_id),
                          l_gene_universe = l_gene_universe,
                          m_go_ids = m_go_ids,
                          l_gene_2_GO = l_gene_2_GO,
                          go_ontology = "BP") %>% 
    mutate(module = .x, .before = 1) 
  
) %>% 
  clean_names() %>%
  dplyr::select(-term) %>%
  left_join(df_go_defs, by = "go_id") %>%
  dplyr::select(module, go_id, term, everything()) %>%
  distinct()

library(RColorBrewer)
library(tidytext)
df_mods_go %>% 
  
  group_by(module) %>% 
  arrange(module, p_adj) %>% 
  dplyr::top_n(n = 3) %>% 
  
  ggplot(aes(x = -log10(p_adj), y = reorder_within(str_wrap(term, width = 25), -p_adj, module))) +
  geom_point(aes(size = significant/annotated, color = -log10(p_adj))) +
  labs(y = "") +
  facet_wrap(~module, scales = "free") +
  scale_y_discrete(labels = function(x) gsub("\\_.*", "", x)) +
  
  geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
  geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
  geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +
  scale_color_gradientn(colors = brewer.pal(9, "Reds")[3:9]) +
  theme(plot.title = element_text(size = 18),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) +
  ggtitle("Gene Ontology enrichment")




