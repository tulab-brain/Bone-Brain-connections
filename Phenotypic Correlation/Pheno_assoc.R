## input:
# - preprocessed bone measure phenotype data file (RDS format): bone_measure_data
# - preprocessed brain measure phenotype data file (RDS format): brain_measure_data
# - covariate information (RDS format): cov_data
## output:
# - regression results for each bone-brain pair

library(dplyr)

pheno_assoc <- function(bone_measure_data, brain_measure_data, cov_data){
  ## load data
  bone2use <- readRDS(brain_measure_data)
  brain2use <- readRDS(bone_measure_data)
  cov2use <- readRDS(cov_data)
  
  ## merge data
  full_data <- merge(cov2use, brain2use, by = "eid") %>% 
    merge(bone2use, by = "eid") %>% 
    as.data.frame()

  ## define variables
  independent_vars <- names(bone2use)[-1]
  dependent_vars <- names(brain2use)[-1]
  covariates <- names(cov2use)[-1]
  
  ## run regression and collect results
  results <- do.call(rbind, lapply(dependent_vars, function(dep_var) {
    do.call(rbind, lapply(independent_vars, function(indep_var) {
      model <- lm(as.formula(paste(dep_var, "~", paste(c(indep_var, covariates), collapse = "+"))),
                  data = full_data)
      result <- as.data.frame(t(summary(model)$coefficients[2, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")] ))
      result$dependent_var <- dep_var
      result$independent_var <- indep_var
      return(result)
    }))
  }))
  
  colnames(results) <- c("Beta", "SE", "t_value", "P_value", "Brain_Measure", "Bone_Measure")
  
  return(results)
}
