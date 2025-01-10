


# libraries ---------------------------------------------------------------


library(tidyverse)
library(tidymodels)
library(tidylo)
library(stm)
library(purrr)
library(RColorBrewer)
library(ggwordcloud)
library(ggrepel)
library(patchwork)

#devtools::install_github("juliasilge/tidytext", force = TRUE) 
library(tidytext)

theme_set(theme_bw())


# data --------------------------------------------------------------------

soft_power <- 3
base_dir <- "~/Documents/work/PhD/projects/WGCNA/sgacc_wgcna_grcca/"
prefix <- "20230228_185samples_19kgenes_vst_qSVA123567_MP_RNAbatch_Race_resids"

load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_GO_RES.RDS")) # df_mods_go


# topic modeling ---------------------------------------------------------


# Unnest tokens - each word is a separate row in a tibble
df_go_tokens <- df_mods_go %>% 
  
  # select the module set we're using (there's many in this tibble across different parameter combos)
  filter(mod_set == "3_40_0.97" & p_adj < 0.05) %>% 
  
  # remove any rows where the GO term is NA
  filter(!is.na(term)) %>% 
  
  # we are running this for each module, so that should be our grouping variable
  group_by(module) %>% 
  
  # split the term descript up so each word is a row in the tibble
  unnest_tokens(word, term) %>% 
  
  # remove stop words (uninteresting words like of, from, and, the, etc.)
  anti_join(get_stopwords())

head(df_go_tokens)

# Create a sparse matrix
go_sparse <- df_go_tokens %>%
  
  # count number of times a word appears in each module
  dplyr::count(module, word) %>%
  
  # convert module from factor to numeric
  mutate(module = as.numeric(as.character(module))) %>% 
  
  # filter for words that appear at least 3 times
  # also, in this case, I've removed words that appear in a lot of path names but are effectively meaningless without context
  filter(n > 3 & !(word %in% c("positive", "negative", "regulation", "response"))) %>%
  cast_sparse(module, word, n)

# Fit topic model using Structural Topic Models (http://www.structuraltopicmodel.com/)
set.seed(0306)
topic_model <- stm(go_sparse, 
                   K = 5, # K is the number of 'topics' (clusters)
                   verbose = FALSE)

summary(topic_model)
# TOPIC DESCRIPTIONS:
# lift = frequency divided by frequency in other topics
# FREX weights words by frequency and exclusivity to the topic

# Ta-da!

# But.... how can we determine the optimal number of topics for our dataset?
df_many_models <- tibble(K = seq(3, 10)) %>% 
  mutate(topic_model = map (
    
    .x = K,
    .f = ~ stm(go_sparse, 
               K = .x, 
               verbose = FALSE)
    
  ))

# Now that we’ve fit all these topic models with different numbers of topics, 
# we can explore how many topics are appropriate/good/“best”

df_diagnostics <- df_many_models %>%
  mutate(exclusivity = map(topic_model, exclusivity),
         semantic_coherence = map(topic_model, semanticCoherence, go_sparse),
         residual = map(topic_model, checkResiduals, go_sparse),
         bound =  map_dbl(topic_model, function(x) max(x$convergence$bound)),
         lfact = map_dbl(topic_model, function(x) lfactorial(x$settings$dim$K)),
         lbound = bound + lfact,
         iterations = map_dbl(topic_model, function(x) length(x$convergence$bound)))

# plot diagnostic results
# we want low residuals, high semantic coherence
df_diagnostics %>%
  transmute(K,
            `Lower bound` = lbound,
            Residuals = map_dbl(residual, "dispersion"),
            `Semantic coherence` = map_dbl(semantic_coherence, mean)) %>%
  gather(Metric, Value, -K) %>%
  ggplot(aes(K, Value, color = Metric)) +
  geom_line(linewidth = 1.5, alpha = 0.7, show.legend = FALSE) +
  facet_wrap(~Metric, scales = "free_y") +
  labs(x = "K (number of topics)",
       y = NULL,
       title = "Model diagnostics by number of topics",
       subtitle = "These diagnostics indicate that a good number of topics would be around 60")


# Semantic coherence is maximized when the most probable words in a given topic frequently co-occur together, 
# and it’s a metric that correlates well with human judgment of topic quality. 
# Having high semantic coherence is relatively easy, though, 
# if you only have a few topics dominated by very common words, 
# so you want to look at both semantic coherence and exclusivity of words to topics. 
# It’s a tradeoff. Read more about semantic coherence in the original paper about it (https://dl.acm.org/citation.cfm?id=2145462.

df_diagnostics %>%
  select(K, exclusivity, semantic_coherence) %>%
  unnest(cols = c(exclusivity, semantic_coherence)) %>%
  mutate(K = as.factor(K)) %>%
  ggplot(aes(semantic_coherence, exclusivity, color = K)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(x = "Semantic coherence",
       y = "Exclusivity",
       title = "Comparing exclusivity and semantic coherence",
       subtitle = "Models with fewer topics have higher semantic coherence for more topics, but lower exclusivity")

# In choosing our final topic model, we should go with a K value that is statistically sound,
# but also biologically meaningful and informative

# plot results ------------------------------------------------------------


topic_model <- stm(go_sparse, 
                   K = 4, # K is the number of 'topics' (clusters)
                   verbose = FALSE)

summary(topic_model)


p_topic <- tidy(topic_model, matrix = "frex") %>%
  group_by(topic) %>%
  slice_head(n = 10) %>%
  mutate(topic = paste0("topic ", topic)) %>% 
  mutate(size = row_number()) %>% 
  
  ggplot(aes(label = term)) +
  geom_text_wordcloud(aes(size = size)) +
  scale_size_continuous(range = c(3, 4)) +
  facet_wrap(vars(topic), nrow = 1) +
  theme_classic()

# GENE LIST - TOPIC PROBABILITIES
df_col <- df_modules_filt %>% 
  dplyr::select(module, color) %>% 
  distinct()
colors <- df_col %>% pull(color)
names(colors) <- df_col %>% pull(module)

group_gamma <- tidy(
  topic_model, 
  matrix = "gamma",
  document_names = rownames(go_sparse)
) %>% 
  mutate(module = factor(document)) %>% 
  dplyr::select(module, topic, gamma) %>% 
  mutate(topic = factor(topic))

p_gamma <- group_gamma %>%
  filter(gamma > 0.01 ) %>% 
  ggplot(aes(x = topic, y = gamma, color = module)) +
  geom_point(aes(alpha = gamma)) +
  geom_text_repel(aes(alpha = gamma, label = module), size = 5, 
                  max.overlaps = 50) +
  scale_color_manual(values = colors) +
  
  # geom_vline(aes(xintercept = 1.5), lty = 2, color = "black") +
  # geom_vline(aes(xintercept = 2.5), lty = 2, color = "black") +
  # geom_vline(aes(xintercept = 3.5), lty = 2, color = "black") +
  # geom_vline(aes(xintercept = 4.5), lty = 2, color = "black") +
  # geom_vline(aes(xintercept = 5.5), lty = 2, color = "black") +
  
  labs(y = expression(gamma), x = "") +
  ggtitle("Module strength of association with each topic") +
  theme(legend.position = "none",
        axis.text.x = element_blank())

(p_gamma / p_topic) + plot_layout(heights = c(2, 1))

