\section{OP 2025}

Notebooks for reproducing all figures and analysis in the Transcriptomic responses to endurance exercise training in rats preprint. 

\subsection*{Getting Started}

All analysis notebooks, saved as .ipynb's in analysis scripts, can be run from Google Colab. Colab links are included in every notebook. 

All saved/processed data used for analysis is streamed to the notebooks from [CaltechData](https://data.caltech.edu/).

\subsection*{Notebooks Directory Contents}
 	1) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Initial_RNA_Analysis.ipynb}{Initial RNA Analysis}
		- Initial clustering and filtering of input matrix
		- Generates Figure 1b
	2) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Label_Correction_Reanalysis.ipynb}{Label Correction and Reanalysis}
		- Rescues mislabeled samples then re-filters/clusters the samples
	3) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/scVI_Batch_Correction.ipynb}{scVI Batch Correction}
		- Generates Supplemental Figures 3a-d
	4) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Rat_Individual_Overview.ipynb}{Rat Individual Overview}
		- Visualizes distribution of samples across -omics and individuals
		- Generates Figure 1a and Supplemental Figure 1
	5) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Linear_Regression_RNA.ipynb}{Linear Regression RNA}
		- Concatenates tissue samples by individual rat then runs linear regression
		- Generates Figure 1d and Supplemental Figures 2 and 6
	6) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Linear_Regression_Physiological.ipynb}{Linear Regression Physiological}
		- Concatenates tissue samples by individual rat then runs linear regression
		- Generates Supplemental Figures 8a, 8b, and 9
	7) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Partial_Correlation_Analysis.ipynb}{Partial Correlation Analysis}
		- Concatenates tissue samples by individual rat then runs linear regression
	8) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Linear_Regression_scVI.ipynb}{Linear Regression scVI}
		- Concatenates tissue samples by individual rat then runs linear regression
		- Generates Supplemental Figure 7
	9) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Rat_Omic_Loading.ipynb}{Rat Omic Loading}
		- Loads -omic data and adds metadata
	10) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Linear_Regression_OMIC.ipynb}{Linear Regression Omic}
		- Concatenates tissue samples by individual rat then runs linear regression
		- Generates Supplemental Figures 4, 10a, and 10b
	11) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/ATAC_METHYL_Analysis.ipynb}{ATAC and METHYL Analysis}
		- Generates PCA plots of ATAC and METHYL data
		- Generates Supplemental Figures 5a, 5b, and 11
	12) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/DEseq.ipynb}{DEseq}
		- Differential gene expression analysis
		- Generates Figures 2a-i
	13) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Virus_QC.ipynb}{Virus QC}
		- Virus quality control to remove known contaminates
	14) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Virus_Correlation_Over_Time.ipynb}{Virus Correlation Over Time}
		- Correlates virus quantity over time
		- Generates Figure 3a
	15) \href{https://github.com/pachterlab/OP_2025/blob/main/analysis_scripts/Virus_DEseq.ipynb}{Virus DEseq}
		- Differential virus expression analysis
		- Generates Figures 3b-e

	