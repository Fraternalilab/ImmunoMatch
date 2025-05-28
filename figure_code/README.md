# code to reproduce figures

This directory contains the code (Python/R) used directly to generate the figures.

|File name|Description|Corresponding figures in manuscript| 
|---------|-----------|-----------------------------------|
|`immunomatch_popscomp.R`| Analysis of antibody structures to compare interface surface area contribution from CDR1/2/3 and FWR in the VH and VL domains|Figure 1c| 
|`comparison_of_diff_models.ipynb`| Comparison of test-set prediction accuracies across different machine learning models|Figure 1d,e,f| 
|`immunomatch_vcab_ltype.R`| Sequence identity comparison of κ- and λ-bearing antibody structures |Figure 2a| 
|`immunomatch_ltype_model_kfold_accuracy.ipynb`| ImmunoMatch-κ and ImmunoMatch-λ accuracy comparison|Figure 2c| 
|`spatial_data_pairing_scores_analysis.ipynb` | Comparison of the 'repair score' from [Engblom et al. <i>Science 2023</i>](https://www.science.org/doi/10.1126/science.adf8486) and the pairing scores calculated from ImmunoMatch |Figure3|
|`immunomatch_cancer.R`| Comparison of pairing scores between healthy Naive, GC, Memory B cells, and in leukaemia/lymphoma samples |Figure 4a,d| 
|`immunomatch_jaffe.R`| Inference of clonotype trees in Jaffe et al. dataset and mapping of pairing scores on the trees|Figure 4b|
|`immunomatch_king_tonsils_scores.R`| Comparison of pairing scores between BCR isotypes and H chain mutational levels in the King et al. Tonsil dataset|Figure 4c|
|`immunomatch_therapeutic_abs.R`| Analysis of pairing scores in therapeutic antibodies|Figure 5|
