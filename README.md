# OP_2025

Notebooks for reproducing all figures and analysis in the Transcriptomic responses to endurance exercise training in rats preprint. 

## Getting Started

All analysis notebooks, saved as .ipynb's in analysis scripts, can be run from Google Colab. Colab links are included in every notebook. 

All saved/processed data used for analysis is streamed to the notebooks from [CaltechData](https://data.caltech.edu/).

## Notebooks Directory Contents

 	1) [Initial RNA Analysis](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Initial_RNA_Analysis.ipynb)
		- Initial clustering and filtering of input matrix
		- Generates Figure 1b
	2) [Label Correction and Reanalysis](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Label_Correction_Reanalysis.ipynb)
		- Rescues mislabeled samples then re-filters/clusters the samples
	3) [scVI Batch Correction](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/scVI_Batch_Correction.ipynb)
		- Generates Supplemental Figures 3a-d
	4) [Rat Individual Overview](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Rat_Individual_Overview.ipynb)
		- Visualizes distribution of samples across -omics and individuals
		- Generates Figure 1a and Supplemental Figure 1
	5) [Linear Regression RNA](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Linear_Regression_RNA.ipynb)
		- Concatenates tissue samples by individual rat then runs linear regression
		- Generates Figure 1d and Supplemental Figures 2 and 6
	6) [Linear Regression Physiological](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Linear_Regression_Physiological.ipynb)
		- Concatenates tissue samples by individual rat then runs linear regression
		- Generates Supplemental Figures 8a, 8b, and 9
	7) [Partial Correlation Analysis](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Partial_Correlation_Analysis.ipynb)
		- Concatenates tissue samples by individual rat then runs linear regression
	8) [Linear Regression scVI](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Linear_Regression_scVI.ipynb)
		- Concatenates tissue samples by individual rat then runs linear regression
		- Generates Supplemental Figure 7
	9) [Rat Omic Loading](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Rat_Omic_Loading.ipynb)
		- Loads -omic data and adds metadata
	10) [Linear Regression Omic](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Linear_Regression_OMIC.ipynb)
		- Concatenates tissue samples by individual rat then runs linear regression
		- Generates Supplemental Figures 4, 10a, and 10b
	11) [ATAC and METHYL Analysis](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/ATAC_METHYL_Analysis.ipynb)
		- Generates PCA plots of ATAC and METHYL data
		- Generates Supplemental Figures 5a, 5b, and 11
	12) [DEseq](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/DEseq.ipynb)
		- Differential gene expression analysis
		- Generates Figures 2a-i
	13) [Virus QC](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Virus_QC.ipynb)
		- Virus quality control to remove known contaminates
	14) [Virus Correlation Over Time](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Virus_Correlation_Over_Time.ipynb)
		- Correlates virus quantity over time
		- Generates Figure 3a
	15) [Virus DEseq](https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Virus_DEseq.ipynb)
		- Differential virus expression analysis
		- Generates Figures 3b-e
	