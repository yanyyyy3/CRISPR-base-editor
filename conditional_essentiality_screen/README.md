This folder contatins codes and files for the *SECOND* screen, *conditional* essentiality screen.
Analysis steps are highly similar to the first screen.

For screen analysis:

1. The sequence data of library screening is available with accession number GSE225335. After merging the pair-end reads, run "conditional_essentiality_screen_read_count.py" to obtain the counts of gRNAs. A summary table of counts from all samples is also available in GSE179913.
2. To obtain the differential abundance (depletion score) of gRNAs, run "conditional_essentiality_screen_differential_abundance.R" . The result files "*_QLFTest.csv" are included here.

For applying machine learning:

1. To optimize machine learning model using automated machine learning tool auto-sklearn, run "conditional_essentiality_screen_autosklearn.py".
2. To evaluate and interprete the optimized model from auto-sklearn, run "conditional_essentiality_screen_sklearn_model_interpretation.py" (optimized model has been implemented in the script).


For each Python script, "-h" shows the detailed description, options and example to run the script.
