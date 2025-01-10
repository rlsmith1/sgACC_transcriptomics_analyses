
########################################################################################

# Illustrate the influence of mu and lambda on regulation of gene and module weights

########################################################################################


df_x_hyperparam %>% 
    filter(str_detect(framework, "lambda0")) %>% 
    group_by(ensembl_gene_id) %>% 
    mutate(mean_weight = mean(weight)) %>% 
    
    ggplot(aes(x = mean_weight)) +
    geom_histogram()

# MU
tibble(
    
    x = 1:1000,
    group = c(rep("a", 375), 
              rep("b", 250), 
              rep("c", 375)),
    
    y1 = c(rnorm(375, mean = 5, sd = 2), 
           rnorm(250, mean = 2, sd = 2), 
           rnorm(375, mean = -4, sd = 2)),
    
    y2 = c(rnorm(375, mean = 3, sd = 1.5), 
           rnorm(250, mean = 0, sd = 1.5), 
           rnorm(375, mean = -2, sd = 1.5)),
    
) %>% 
    pivot_longer(starts_with("y"), names_to = "facet", values_to = "y") %>% 
    
    ggplot(aes(x = x, y = y)) +
    geom_col(aes(alpha = y, color = group)) +
    scale_color_manual(values = c("#ffda80ff", "#80ffd7ff", "#fa86f8f5")) +
    
    facet_wrap(vars(facet)) +
    labs(x = "", y = "weights") +
    #ggtitle("shrinkage regularization") +
    theme(strip.text = element_blank(),
          axis.text = element_blank(),
          legend.position = "none")

# LAMBDA
tibble(
    
    x = 1:1000,
    group = c(rep("a", 375), 
              rep("b", 250), 
              rep("c", 375)),
    
    y1 = c(rnorm(375, mean = 5, sd = 2), 
           rnorm(250, mean = 2, sd = 2), 
           rnorm(375, mean = -4, sd = 2)),
    
    y2 = c(rnorm(375, mean = 5, sd = 0.5), 
           rnorm(250, mean = 2, sd = 0.5), 
           rnorm(375, mean = -4, sd = 0.5)),
    
) %>% 
    pivot_longer(starts_with("y"), names_to = "facet", values_to = "y") %>% 
    
    ggplot(aes(x = x, y = y)) +
    geom_col(aes(alpha = y, color = group)) +
    scale_color_manual(values = c("#ffda80ff", "#80ffd7ff", "#fa86f8f5")) +
    
    facet_wrap(vars(facet)) +
    labs(x = "", y = "weights") +
    #ggtitle("shrinkage regularization") +
    theme(strip.text = element_blank(),
          axis.text = element_blank(),
          legend.position = "none")



