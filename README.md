# A neuro-immune axis of transcriptomic dysregulation within the subgenual anterior cingulate cortex in schizophrenia
### Author: Rachel L Smith

This repository contains the code used to run data preprocessing and analysis and generate figures for R.L. Smith et al 2025 (bioRxiv: ).

For details on the CCA/PLS toolkit, refer to the README available in RCCA_toolkit/cca_pls_toolkit_final.

All functions can be found in functions/ and are loaded in scripts/setup.R.

-------------------------------------------------------------------
The scripts directories are organized as follows:
* **00.preprocessing** - Clean/organize study covariates and preprocess gene & transcript raw counts; run MCA on toxicology data
* **01.DGE** - Run standard differential gene expression analysis using DESeq2; characterize results enrichment using GSEA
* **02.WGCNA** - Run weighted gene co-expression network analysis (WGCNA); characterize module enrichments; run standard module eigengene analysis to identify module-covariate associations
* **03.CCA_toolkit_setup** - Format and export data to run GRCCA on the cluster (example code to run the analysis can be found in this folder as well as ~/RCCA_toolkit/cca_pls_toolkit_final)
* **04.GRCCA** - Read GRCCA results (generated on cluster using grcca_analysis.m); characterize enrichment of results using GSEA
* **05.GRCCA_transcripts** - Run WTCNA and GRCCA on transcript-level data (both on cluster); characterize enrichment of results using GSEA
* **06.sensitivity_RCCA** - Read RCCA results (generated on cluster using rcca_analysis.m); characterize results enrichment using GSEA and compare to GRCCA
* **07.figures_markdown** - R markdown files to generate figures 1-5, as well as all supplementary figures
* Other scripts (not organized into subfolder)
    - **load_genes.R** - Loads preprocessed gene-level data (run after 00.preprocessing but before 02.WGCNA)
    - **load_transcripts.R** - Loads preprocessed transcript-level data and gene ID mappings (run after 00.preprocessing but before 05.GRCCA_transcripts)
    - **load_WGCNA_res.R** - Loads WGCNA module assignments, as well as data from load_genes.R (run after 02.WGCNA but before 03.CCA_toolkit_setup)
    - **load_WTCNA_res.R** - Loads WTCNA module assignments, as well as data from load_transcripts.R (run after 05.GRCCA_transcripts/00.case_control_WTCNA.R)
    - **setup.R** - Loads all necessary packages, study covariate data, cell type data, gene ontology pathways, risk gene lists; sources functions in ~/functions; sets plot theme for figures in 07.figures_markdown (run before running anything else!)
-------------------------------------------------------------------

Individual scripts are commented to describe their specific purpose. For any questions on code/implementation, please email corresponding author rachel.smith2@nih.gov.
