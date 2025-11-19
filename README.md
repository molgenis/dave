# Digital Approximation of Variant Effect (MOLGENIS DAVE)

Genome diagnostics remain limited by the high fraction of variants of uncertain significance (VUS), largely due to insufficient interpretability of missense variation. Existing pathogenicity predictors offer strong performance but often lack transparency or mechanistic insight. Here, we present the Digital Approximation of Variant Effect, version 1 (DAVE), an explainable missense variant predictor built on 12 biophysically grounded features spanning stability, hydrophobicity, electrostatics, and molecular interactions. Trained on curated Dutch diagnostic data, DAVE achieves reliable classification while decomposing predictions into interpretable feature contributions. This framework enables mechanistic follow-up, reduces VUS burden, and advances clinically actionable variant interpretation.

### Publication
* _in writing_

### Key resources
* [Freeze 6 dataset](data/freeze6.csv.gz) (Gzipped CSV file)
  * 23417 rows, each a DNA/protein variant (9345 LB, 2703 LP, 11221 VUS and 148 CF)
  * 81 columns, metadata and features
* [Modeling and prediction on freeze 6](src/main/R/ML_on_freeze6.R) (R script)
* [Random Forest prediction model](models/dave_rf_model.RData) (RData file)

### Raw computational results
* [Wild-type and mutant proteins](data/genes)
