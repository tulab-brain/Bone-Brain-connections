library(data.table)
library(dplyr)
library(stringr)
library(TwoSampleMR)
library(RadialMR)
library(mr.raps)
library(ieugwasr)
library(plinkbinr)

####Define Function####
source("/mr_modified.R") # The mr-raps cannot not be conducted in the original function; derived from https://github.com/linjf15/MR_tricks 
source("/MRFindEAF.R") # Find effect allele frequency; derived from https://github.com/linjf15/MR_tricks 

Prepro <- function(data, threshold_GWAS) {
  if ("pval.exposure" %in% names(data)) {
    data <- data[data$MAF >= 0.01 & (data$beta.exposure^2 / data$se.exposure^2) > 10, ] # Reduce weak instrument bias
  } else if ("pval.outcome" %in% names(data)) {
    data <- data[data$pval.outcome > threshold_GWAS, ]  # remove the variants strongly associated with outcomes, reducing false positive
  }
  return(data)
} 

Plink <- get_plink_exe() # Get LD matrix using local plink binary and reference dataset

####Define path####
Exposure_Dir <- "path to Exposure files" 
Outcome_Dir <- "path to Outcome files" 
Output_Dir <- 'path to save MR estimate results' 
ScalePara <-  readRDS('Path to Scale files') # Scale beta and SE of bone measures based on phenotypes'SD

#### Define Parameter####
threshold_GWAS <- 5e-8  # Threshold for IVs

####Run####
filenamesBone <- list.files(path = Exposure_Dir, pattern = "*.txt")
filenamesBrainTraits <- list.files(path = Outcome_Dir, pattern = "*.txt")

epoch = 1
for (exposure_index in filenamesBone)  {
  exposure <- fread(str_c(Exposure_Dir,exposure_index))
  setnames(exposure, new = c('pval'), old = c('P_value'))
  exposure <- exposure[which(exposure$pval < threshold_GWAS),]
  exposure$rsid <- exposure$ID 
  
  ## Scale Beta and SE of bone measures
  Scale_SD <- subset(ScalePara,ScalePara==str_extract(exposure_index,"\\d+"))$SD
  exposure$BETA <- exposure$BETA/Scale_SD; exposure$SE = exposure$SE/Scale_SD;
  
  ## Clumping
  exposure <- ld_clump(
    exposure,
    clump_kb = 1000,
    clump_r2 = 0.001,
    plink_bin = Plink,
    bfile = "/1kgClump/EUR" #use local file in case of network block
    ) 
  
  exposure <- as.data.frame(exposure)
  
  exposure <- format_data(
    dat= exposure,
    type= "exposure",
    header = TRUE,
    snp_col = "rsid",
    beta_col =  "BETA",
    se_col = "SE",
    pval_col = "pval",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    chr_col = "CHROM",
    pos_col = "GENPOS"
  )

  exposure$MAF <- ifelse(exposure$eaf.exposure > .5, 1-exposure$eaf.exposure, exposure$eaf.exposure)
  exposure <- Prepro(exposure,threshold_GWAS)
  
  for (outcome_index in filenamesBrainTraits)  {
    Outcome <- fread(str_c(Outcome_Dir,outcome_index),fill = TRUE)
    Outcome <- Outcome[Outcome$snpid %in% exposure$SNP,]
    Outcome$pval <- as.numeric(Outcome$pval)
    Outcome <- as.data.frame(Outcome)
    Outcome <- format_data(
      dat= Outcome,
      type= "outcome",
      header = TRUE,
      snp_col = "snpid",
      beta_col =  "beta",
      se_col = "se",
      pval_col = "pval",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      eaf_col = "A1fre"
    )
    
    if (any(is.na(Outcome$eaf.outcome))){Outcome = snp_add_eaf(Outcome)}

    Outcome <- Prepro(Outcome,threshold_GWAS)
    
    ## Harmonise
    mydata <- harmonise_data(
      exposure_dat= exposure,
      outcome_dat= Outcome,
      action= 2)
  
    ## Remove outliers using RadialMR
    outliers <- ivw_radial(r_input = mydata, alpha = 0.05, weights = 1, tol = 0.0001, summary = TRUE)
    mydata <- mydata[!(mydata$SNP %in% outliers[["outliers"]][["SNP"]]),]

    outliers <- egger_radial(r_input = mydata, alpha = 0.05, weights = 1, summary = TRUE)
    mydata <- mydata[!(mydata$SNP %in% outliers[["outliers"]][["SNP"]]),] 
    
    
    ##  IVW-random effect when IVs more than 3, otherwise -fix effect
    if (nrow(mydata)>3){ res <- mr_modified(mydata,method_list = c("mr_ivw_mre","mr_egger_regression", "mr_raps",
                                                                   "mr_simple_median","mr_weighted_median","mr_penalised_weighted_median",
                                                                   "mr_simple_mode","mr_weighted_mode","mr_simple_mode_nome","mr_weighted_mode_nome") )
    }    else if (nrow(mydata)>1) {

   
      res <- mr_modified(mydata,method_list = c("mr_ivw_fe","mr_egger_regression", "mr_raps",
                                                "mr_simple_median","mr_weighted_median","mr_penalised_weighted_median",
                                                "mr_simple_mode","mr_weighted_mode","mr_simple_mode_nome","mr_weighted_mode_nome"))
    } else {
  
      res <- mr_modified(mydata,method_list = c("mr_wald_ratio"))} # note: MR estimate with only 1 IV will not be considered
    
    epoch <- epoch+1
  }
}