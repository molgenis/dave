# Digital Approximation of Variant Effect (MOLGENIS DAVE)

Diagnostic yield in NGS genome diagnostics is constraint by the high fraction of variants of uncertain significance (VUS), in large part due to insufficient interpretability of missense variation. Existing pathogenicity predictors offer strong performance, but often produce an unexplainable score lacking mechanistic insight. Here, we present the Digital Approximation of Variant Effects (MOLGENIS DAVE), an explainable missense variant predictor built on 12 biophysically grounded features spanning stability, hydrophobicity, electrostatics, and molecular interactions. Trained on curated Dutch diagnostic data, DAVE reliably classifies and breaks down predictions into interpretable feature contributions. With a focus on explainability, this framework aims to alleviate the VUS burden, advances clinically actionable variant interpretation and enables mechanistic follow-up.

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
