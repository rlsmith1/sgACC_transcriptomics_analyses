
###############################################################################################################

## Access the most recent release of the Ensembl database and convert to a tibble for easy use in other scripts

###############################################################################################################

## As of 11-Dec-2023, most current release of GRCh38 is version 110 (2023-05-04)
## Downloaded from https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/

## LIBRARIES
library(tidyverse)
library(rtracklayer)
library(Repitools)

## DATA
#hsapiens_path <- "~/Downloads/Homo_sapiens.GRCh38.110.chr.gtf"
hsapiens_path <- "~/Downloads/Homo_sapiens.GRCh38.110.gtf"
hsapiens_data <- rtracklayer::import(hsapiens_path)

# Convert to dataframe and save for future scripts
df_hsapiens_genome <- annoGR2DF(hsapiens_data) %>% 
  as_tibble() %>% 
  mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name),
         gene_name = ifelse(is.na(gene_name), gene_id, gene_name)
  )

## SAVE 
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
save(df_hsapiens_genome,
     file = paste0(base_dir, "objects/hsapiens_genome_v110.RDS")
)

## ALTERNATIVE
# can also access these data using biomaRt:

# library(biomaRt)
# hsapiens_mart <- useEnsembl(biomart = "ensembl", 
#                                dataset = "hsapiens_gene_ensembl",
#                                 version = 110
# )



