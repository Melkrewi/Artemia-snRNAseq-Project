### Analysis scripts:
run_analysis_Harmony_DUBStepR.R: This Rscript takes as input the output of cellbender for the different replicates and the filtered cells (from QC), and performs the clustering analysis.
atac_fragments

specificity_analysis.R: This script performs the marker specificty analysis used to choose the appropriate clustering resolution. 

new_atac_analysis_v3.R: This performs the analysis and integration of the ATACseq data from the multiome-replicates. Peaks were called for replicate 3 using callpeaks_replicate_3_v3.R and for replicate 4 using callpeaks_replicate_4_v3.R.

pseudotime_multiome.R: This script performs the pseudotime analysis on the output of new_atac_analysis_v3.R. 



