# A transcriptomic dimension of neuronal and immune gene programs within the subgenual anterior cingulate cortex in schizophrenia
### Author: Rachel L Smith

This repository contains the code used to run data preprocessing and analysis and generate figures for R.L. Smith et al Transl Psych (2026).

For details on the CCA/PLS toolkit, refer to the README available in RCCA_toolkit/cca_pls_toolkit_final.

All functions can be found in functions/ and are loaded in code/analysis/setup.R.

-------------------------------------------------------------------
The scripts directories are organized as follows:
* **00_preprocessing** - Clean/organize study covariates and preprocess gene & transcript raw counts; run MCA on toxicology data
* **01_WGCNA** - Run weighted gene co-expression network analysis (WGCNA); characterize module enrichments; run standard module eigengene analysis to identify module-covariate associations
* **02a_CCA-toolkit-setup** - Format and export data to run GRCCA on the cluster (example code to run the analysis can be found in this folder as well as ~/RCCA_toolkit/cca_pls_toolkit_final)
* **02b_GRCCA-genes** - Read GRCCA results (generated on cluster using grcca_analysis.m); characterize enrichment of results using GSEA; sensitivity analyses
* **02c_GRCCA-transcripts** - Run WTCNA and GRCCA on transcript-level data (both on cluster); characterize enrichment of results using GSEA
* **03_DGE** - Run standard differential gene expression analysis using limma-voom; characterize results enrichment using GSEA
* Helper scripts (not organized into subfolder)
    - **load_genes.R** - Loads preprocessed gene-level data (run after 00_preprocessing but before 02_WGCNA)
    - **load_transcripts.R** - Loads preprocessed transcript-level data and gene ID mappings (run after 00_preprocessing but before 03d_GRCCA-transcripts)
    - **load_WGCNA_res.R** - Loads WGCNA module assignments, as well as data from load_genes.R (run after 01_WGCNA but before 02a_CCA-toolkit-setup)
    - **load_WTCNA_res.R** - Loads WTCNA module assignments, as well as data from load_transcripts.R (run after 02c_GRCCA-transcripts/00_case-control-WTCNA.R)
    - **setup.R** - Loads all necessary packages, study covariate data, cell type data, gene ontology pathways, risk gene lists; sources functions in ~/functions; sets plot theme for figures
-------------------------------------------------------------------

Individual scripts are commented to describe their specific purpose. For any questions on code/implementation, please email corresponding author at smith.rachel.lillian@gmail.com.


![Graphical abstract](outputs/figures/main/00_combined_panels/Fig1.jpg)
