
################################################################################################################

# TRANSCRIPTS: Filter transcripts based on coefficient of variation; i.e., how variable is the transcript across
# samples given its expression levels?

################################################################################################################


# Setup -------------------------------------------------------------------

## Set directories
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
supplement_figures_dir <- paste0(base_dir, "outputs/figures/supplement/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

## Load regressed transcript data (all transcripts)
prefix <- "08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC"
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress (generated in 03.transcript_raw_count_preprocessing.R)

## Load raw transcript counts
df_transcript_raw_counts <- read.csv(paste0(base_dir, "data/transcript_count_matrix.csv")) %>% 
    as_tibble() %>% 
    dplyr::rename("ensembl_transcript_id" = "X") %>% 
    rename_all(~str_remove(.x, "_MergedBam_stringtieOutput")) %>% 
    clean_names() %>% 
    rename_all(~str_remove(.x, "x"))

# Calculate transcript coefficient of variation ----------------------------------------------


## "Selecting information features and reducing dimensionality"
# https://www.nature.com/articles/s41576-023-00586-w
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-193 (pick arbitrary threshold?)

# However, we don't just want to filter on variance because variance scales with mean and some transcripts have very low expression
# So we can use coef of variance, which filters on variance relative to mean


### Filter by sd/mean (CV) of raw counts across samples ###

## Calculate CV of each transcript
df_mean_sd <- df_transcript_raw_counts %>% 
    filter(ensembl_transcript_id %in% colnames(df_vsd_regress)) %>% 
    pivot_longer(2:ncol(.), names_to = "sample", values_to = "raw_counts") %>% 
    group_by(ensembl_transcript_id) %>% 
    summarise(
        mean = mean(raw_counts, na.rm = TRUE),
        stdev = sd(raw_counts, na.rm = TRUE),
        CV = stdev/mean
    ) %>% 
    arrange(-CV)
df_mean_sd <- df_mean_sd %>% mutate(z_CV = scale(CV)[,1])

## Determine the CV cutoff based on transcript data spread
# low CV means the transcript expression doesn't really change across samples, thus it will not be informative in our GRCCA analysis (or WGCNA for that matter)
quantiles <- summary(df_mean_sd %>% pull(CV))
cutoff <- quantiles[2]


# Filter transcripts based on desired CV cutoff ------------------------------------------------------


## Keep transcripts that are above CV cutoff
keep <- df_mean_sd %>% filter(CV >= cutoff) %>% pull(ensembl_transcript_id)
length(keep) # n = 54302 (CVq1; where q1 ~= 0.36)

df_vsd_regress_filt <- df_vsd_regress %>% dplyr::select(sample, any_of(keep)) 
ncol(df_vsd_regress_filt) - 1 # n = 54302



# Save objects for downstream analyses ------------------------------------

## Save CV objects for supplemental figures
save(
    df_mean_sd, cutoff,
    file = paste0(analysis_objects_dir, "transcript_coefficient_of_variation.RData")
)


## Save df_vsd_regress_filt for downstream analyses
prefix2 <- paste0(prefix, "_CVq1")
save(df_vsd_regress_filt, file = paste0(base_dir, "objects/", prefix2, ".RDS")) # df_vsd_regress_filt

## Save filtered and regressed counts as a matrix to run WTCNA on the cluster
count_matrix <- df_vsd_regress_filt %>% column_to_rownames("sample")
save(count_matrix, file = paste0(base_dir, "objects/", prefix2, "_matrix.RDS"))


