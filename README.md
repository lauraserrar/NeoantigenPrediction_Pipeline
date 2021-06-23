# NeoantigenPrediction_Pipeline

## Description
This pipeline intends to provide an integrative computational pipeline for neoantigen prediction. Offering the possibility to integrate mutation calling, gene expression, HLA typing and neoantigen prediction in a single computational tool condensing multiple pipelines to provide a final list of significant neoantigens starting only from the raw genome sequencing data files. 

This repository contains most of the scripts used during my final degree project (Apr-Jun 2021)

## Contents

- Scripts:

Contains the scripts for predicting neoantigens.

- Results:
Contains all the outputs of the pipeline per patient.
	- Final folder contains the Additional Supplementary Tables with the final outputs of the pipeline.

## Requirements

In order to run the pipeline additional files are needed.
- List of patients, must be saved in the general folder
- List of samples (per each patient), must be saved in the patient folder
- List of tumour sample (per each patient)must be saved in the patient folder
- List of normal sample (per each patient)must be saved in the patient folder
- List of chromosomes (be careful if X or Y needed) must in saved at the general folder

## About this repository

We summarize the content of the repository here:

```
NeoantigenPrediction_Pipeline
│   README.md
│
└───Scripts
│      run_pyclone.py
│      pipeline.sh
│      sequenza_results.R
│
└───Results
│   │
│   └───Patient_288
│	│   
│   	└───NeoPrePipe
│		288_005_somatic_filtered_PASS.vcf
│		288_005.neoantigens.Indels.summarytable.txt
│		288_005.neoantigens.Indels.txt
│		288_005.neoantigens.summarytable.txt
│		288_005.neoantigens.txt
│		288_006_somatic_filtered_PASS.vcf
│		288_006.neoantigens.Indels.summarytable.txt
│		288_006.neoantigens.Indels.txt
│		288_006.neoantigens.summarytable.txt
│		288_006.neoantigens.txt
│		288_008_somatic_filtered_PASS.vcf
│		288_008.neoantigens.Indels.summarytable.txt
│		288_008.neoantigens.Indels.txt
│		288_008.neoantigens.summarytable.txt
│		288_008.neoantigens.txt
│		288_0016_somatic_filtered_PASS.vcf
│		288_0016.neoantigens.Indels.summarytable.txt
│		288_0016.neoantigens.Indels.txt
│		288_0016.neoantigens.summarytable.txt
│		288_0016.neoantigens.txt
│	│   
│   	└───Optitype
│		2021_06_16_08_55_45_coverage_plot.pdf
│		2021_06_16_08_55_45_result.tsv
│	│   
│   	└───Quantiseq
│		288-005_cell_fractions.txt
│		288-005_gene_count.txt
│		288-005_gene_tpm.txt
│		288-006_cell_fractions.txt
│		288-006_gene_count.txt
│		288-006_gene_tpm.txt
│		288-008_cell_fractions.txt
│		288-008_gene_count.txt
│		288-008_gene_tpm.txt
│		288-016_cell_fractions.txt
│		288-016_gene_count.txt
│		288-016_gene_tpm.txt
│	│   
│   	└───Sequenza
│	    │   
│   	    └───288_005
│		   288_005_alternative_fit.pdf
│		   288_005_alternative_solutions.txt
│		   288_005_chromosome_depths.pdf
│		   288_005_chromosome_view.pdf
│		   288_005_CN_bars.pdf
│		   288_005_confints_CP.txt
│		   288_005_CP_contours.pdf
│		   288_005_gc_plots.pdf
│		   288_005_genome_view.pdf
│		   288_005_model_fit.pdf
│		   288_005_alternative_solutions.txt
│		   288_005_mutations.txt
│		   288_005_segments.txt
│		   288_005_sequenza_cp_table.RData
│		   288_005_sequenza_extract.RData
│		   288_005_sequenza_log.txt
│	    │   
│   	    └───288_006
│		   288_006_alternative_fit.pdf
│		   288_006_alternative_solutions.txt
│		   288_006_chromosome_depths.pdf
│		   288_006_chromosome_view.pdf
│		   288_006_CN_bars.pdf
│		   288_006_confints_CP.txt
│		   288_006_CP_contours.pdf
│		   288_006_gc_plots.pdf
│		   288_006_genome_view.pdf
│		   288_006_model_fit.pdf
│		   288_006_alternative_solutions.txt
│		   288_006_mutations.txt
│		   288_006_segments.txt
│		   288_006_sequenza_cp_table.RData
│		   288_006_sequenza_extract.RData
│		   288_006_sequenza_log.txt
│	    │   
│   	    └───288_008
│		   288_008_alternative_fit.pdf
│		   288_008_alternative_solutions.txt
│		   288_008_chromosome_depths.pdf
│		   288_008_chromosome_view.pdf
│		   288_008_CN_bars.pdf
│		   288_008_confints_CP.txt
│		   288_008_CP_contours.pdf
│		   288_008_gc_plots.pdf
│		   288_008_genome_view.pdf
│		   288_008_model_fit.pdf
│		   288_008_alternative_solutions.txt
│		   288_008_mutations.txt
│		   288_008_segments.txt
│		   288_008_sequenza_cp_table.RData
│		   288_008_sequenza_extract.RData
│		   288_008_sequenza_log.txt
│	    │   
│   	    └───288_016
│		   288_016_alternative_fit.pdf
│		   288_016_alternative_solutions.txt
│		   288_016_chromosome_depths.pdf
│		   288_016_chromosome_view.pdf
│		   288_016_CN_bars.pdf
│		   288_016_confints_CP.txt
│		   288_016_CP_contours.pdf
│		   288_016_gc_plots.pdf
│		   288_016_genome_view.pdf
│		   288_016_model_fit.pdf
│		   288_016_alternative_solutions.txt
│		   288_016_mutations.txt
│		   288_016_segments.txt
│		   288_016_sequenza_cp_table.RData
│		   288_016_sequenza_extract.RData
│		   288_016_sequenza_log.txt
│	│   
│   	└───Final
│	    Final_288_005
│	    Final_288_006
│	    Final_288_008
│	    Final_288_016
│   │   
│   └───Patient_323
│	│   
│   	└───Optitype
│		2021_06_14_11_36_22_result.tsv
│		2021_06_14_11_36_22_coverage_plot.pdf
│	│   
│   	└───Sequenza
│	    │   
│   	    └───323_003
│		   323_003_alternative_fit.pdf
│		   323_003_alternative_solutions.txt
│		   323_003_chromosome_depths.pdf
│		   323_003_chromosome_view.pdf
│		   323_003_CN_bars.pdf
│		   323_003_confints_CP.txt
│		   323_003_CP_contours.pdf
│		   323_003_gc_plots.pdf
│		   323_003_genome_view.pdf
│		   323_003_model_fit.pdf
│		   323_003_alternative_solutions.txt
│		   323_003_mutations.txt
│		   323_003_segments.txt
│		   323_003_sequenza_cp_table.RData
│		   323_003_sequenza_extract.RData
│		   323_003_sequenza_log.txt
│	    │   
│   	    └───323_004
│		   323_004_alternative_fit.pdf
│		   323_004_alternative_solutions.txt
│		   323_004_chromosome_depths.pdf
│		   323_004_chromosome_view.pdf
│		   323_004_CN_bars.pdf
│		   323_004_confints_CP.txt
│		   323_004_CP_contours.pdf
│		   323_004_gc_plots.pdf
│		   323_004_genome_view.pdf
│		   323_004_model_fit.pdf
│		   323_004_alternative_solutions.txt
│		   323_004_mutations.txt
│		   323_004_segments.txt
│		   323_004_sequenza_cp_table.RData
│		   323_004_sequenza_extract.RData
│		   323_004_sequenza_log.txt
```







