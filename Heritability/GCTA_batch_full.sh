#!/bin/bash
 pheno_folder="/media/lei/lei/BoneBrain/Block2_Genetic_Basis/Main/Step1_Heritability/Data/FIRST_GCTA"
 covariate_path="/media/lei/lei/BoneBrain/Block2_Genetic_Basis/Main/Step1_Heritability/Data/cov_bone_gwas.txt"
 GRM_path="/media/lei/lei/BoneBrain/Block2_Genetic_Basis/Main/Step1_Heritability/Data/GCTA_mtx_full/qc_bone_imp"
 output_folder="/media/lei/lei/BoneBrain/Block2_Genetic_Basis/Main/Step1_Heritability/Results Summary/GCTA/FIRST_full"
 
 for pheno in "${pheno_folder}"/*.txt; do
    pheno_name=$(basename "{$pheno}" .txt)
    output_prefix="${output_folder}/GREML_${pheno_name}"
    ./gcta-1.94.1 \
      --reml \
      --grm "$GRM_path" \
      --pheno "${pheno}" \
      --qcovar "${covariate_path}" \
      --out "$output_prefix" \
      --thread-num 85
 done
