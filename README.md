# OP_2025

Notebooks for reproducing all figures and analysis in the Transcriptomic responses to endurance exercise training in rats preprint. 

## Getting Started

All analysis notebooks, saved as .ipynb's in analysis scripts, can be run from Google Colab. Colab links are included in every notebook. 

All saved/processed data used for analysis is streamed to the notebooks from [CaltechData](https://data.caltech.edu/).

## Notebooks Directory Contents

1. [Initial RNA Analysis](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/1_Initial_RNA_Analysis.ipynb)
	- Loads raw RNA counts (879 cells × 24,770 genes) and merges with phenotype metadata
	- Flags UMI outliers using a ±3 SD threshold per tissue; normalizes, log-transforms, and selects highly variable genes
	- Runs PCA and Leiden clustering; assigns tissue labels by majority vote within each cluster and flags candidate mislabelings
	- Generates Figure 1b

2. [Label Correction and Reanalysis](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/2_Label_Correction_Reanalysis.ipynb)
	- Manually corrects two mislabeled samples identified via PCA/cluster mismatch in notebook 1 and removes three highly suspicious samples
	- Propagates corrected metadata to transcript-level data (72,032 features) via barcode matching

3. [scVI Batch Correction](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/3_scVI_Batch_Correction.ipynb)
	- Trains a scVI variational autoencoder (tissue as batch key, sex as covariate, 400 epochs) on normalized counts
	- Extracts a 10-dimensional latent representation and visualizes via PCA; batch correction success assessed visually
	- Generates Supplemental Figure 3

4. [Rat Individual Overview](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/4_Rat_Individual_Overview.ipynb)
	- Loads all 9 omic datasets and builds a binary tissue × individual coverage matrix for each omic
	- Visualizes sample coverage as heatmaps by individual (50 rats) and by tissue (18 tissues × 2 sexes)
	- Generates Figure 1a and Supplemental Figures 1 and 2

5. [Linear Regression RNA](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/5_Linear_Regression_RNA.ipynb)
	- Reshapes RNA data into an individuals-as-rows matrix by concatenating genes across tissues (~286k features); library-size normalizes without log-transforming to satisfy linear regression assumptions
	- Fits Ridge regression to predict weeks of exercise using a ~1/3 train split stratified by time × sex; evaluates with concordance correlation coefficient (CCC) across all tissues combined and each tissue individually
	- Runs pathway enrichment on the top model-weight genes (positive and negative) via rat-to-human gene mapping and GSEA
	- Generates Figure 1d and Supplemental Figures 5, 6, 13, 14–30, and 31

6. [Linear Regression Physiological](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/6_Linear_Regression_Physiological.ipynb)
	- Uses the same reshape and split pipeline as notebook 5 to predict physiological outcomes: % body fat change and VO₂max change
	- Fits separate Ridge regression models for each trait; available sample sizes are smaller (20–23 individuals) due to missing phenotype measurements
	- Generates Supplemental Figures 8 and 9

7. [Partial Correlation Analysis](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/7_Partial_Correlation_Analysis.ipynb)
	- Takes the top-ranked gene from notebook 5 and computes a 4×4 partial correlation matrix with weeks of exercise, % body fat change, and VO₂max change across 29 individuals using `pingouin.partial_corr`
	- Results are exploratory: gene selection is post-hoc based on regression coefficient magnitude

8. [Linear Regression scVI](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/8_Linear_Regression_scVI.ipynb)
	- Reshapes RNA data identically to notebook 5, trains scVI on the reshaped matrix, then fits Ridge regression on the 10D latent space rather than raw gene features to predict weeks of exercise
	- Generates Supplemental Figure 4

9. [Rat Omic Loading](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/9_Rat_Omic_Loading.ipynb)
	- Loads 7 raw omic h5ad files (ATAC, PROT, PHOSPHO, UBIQ, METAB, IMMUNO, ACETYL) and annotates each with standardized metadata and gene/feature identifiers
	- ATAC data receives ChipSeeker annotations (promoter, intron, distance to TSS); methylation requires merging three partial files (~1.59M features total); missing values are dropped rather than imputed

10. [Linear Regression Omic](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/10_Linear_Regression_OMIC.ipynb)
	- Concatenates all 9 omics for Heart tissue only (the only tissue sampled across all omics) into a multi-omic feature matrix; restricts to female individuals due to sparse male coverage at some timepoints
	- Fits Ridge regression per-omic and on the combined matrix to predict weeks of exercise; compares CCC scores across omics
	- Generates Supplemental Figures 10 and 32–40

11. [ATAC and METHYL Analysis](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/11_ATAC_METHYL_Analysis.ipynb)
	- Filters ATAC to promoter peaks, matches to shared genes with RNA, and computes Pearson correlation across 12.3M gene-individual-tissue pairs (r=0.108)
	- Runs PCA on the full methylation feature matrix (~1M sites) and visualizes by tissue, sex, and time
	- Generates Supplemental Figures 7 and 11

12. [DEseq](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/12_DEseq.ipynb)
	- Runs DESeq2 on gene and transcript counts with design formula `~(tissue + sex) * time`, capturing tissue- and sex-specific time responses; features with fewer than 10 total counts are filtered
	- Main contrast: SKM-GN week 8 vs. week 0; results merged with gene symbols and GO annotations
	- Generates Figures 2a-i

13. [Virus Filtering](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/13_Virus_Filtering.ipynb)
	- Loads 99,228 viral features from palmdb quantification; removes known lab contaminants (viruses present in negative control samples) and applies an abundance filter (≥200 cells), retaining 1,492 viruses
	- Assigns ICTV taxonomy and generates interactive Krona plots; ~30k viruses with no taxonomic assignment are retained for downstream analysis

14. [Virus Correlation Over Time](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/14_Virus_Correlation_Over_Time.ipynb)
	- Aggregates total virus counts per individual across 13 shared tissues (39 individuals) and computes Pearson correlation with weeks of exercise (r=0.24, NS) and animal lifetime (r=0.49, p=0.001)
	- Lifetime encodes time=0 controls as 8 weeks, placing all animals on a common survival timeline rather than an exercise-duration timeline
	- Generates Figure 3a

15. [Virus DEseq](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/15_Virus_DEseq.ipynb)
	- Mirrors notebook 12 on the 1,492 filtered viruses using the same `~(tissue + sex) * time` design and week 8 vs. week 0 contrast; tissues causing convergence issues (VENACV, TESTES, OVARY) are excluded
	- Generates Figures 3b-e

16. [Virus BLAST](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/16_Virus_BLAST.ipynb)
	- Extracts raw reads for a virus of interest, runs `blastn` against NCBI nt (top 20 hits per read), removes polyA/T sequencing artifacts (≥12 consecutive A/T bases), and maps hits to NCBI taxonomy
	- Results are illustrative: only 10 reads are sampled per virus, sufficient for identity confirmation but not abundance estimation
	- Generates Supplemental Figure 12
	