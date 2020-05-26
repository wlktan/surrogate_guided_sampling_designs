# surrogate_guided_sampling_designs
This repository contains code for the paper 'Surrogate-guided sampling designs for classification of rare outcomes from electronic medical records data'. The original dataset used in the paper is not available for public distribution, so we include a derived dataset (created text features, instead of original radiology report text) that can be used to reproduce the analysis.

R file descriptions:
- `display_results.R`: Displays results as in Table 1 of the manuscript from pre-loaded data analysis results (i.e. `p1_lasso_n100.csv`, `p1_lasso_n250.csv`, `p1_lasso_n500.csv` files)

- `run_data_example_analysis.R`: Runs the data analysis as described in Section 4 of the manuscript, using derived dataset `sgs_data_application_feature_matrix.csv`)

CSV file descriptions:
- `p1_lasso_n100.csv`, `p1_lasso_n250.csv`, `p1_lasso_n500.csv`: The output from `run_data_example_analysis.R` - validation AUC results for training sample sizes 100, 250, and 500, and sampling using SRS and SGS (columns) over 1000 iterations (one per row).

- `sgs_data_application_feature_matrix.csv`: Derived dataset from the original annotated LIRE dataset, where the unigram (single word) features were created with the term-frequency inverse-document frequency (TF-IDF) representation, excluding typical English stopwords as well as terms that were very rare (<5% of all reports) or common (>90% of all reports).

