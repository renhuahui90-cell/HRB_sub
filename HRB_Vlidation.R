#------------------loading package-------------------------------# 
setwd("./")
rm(list=ls())

#library(rmeta)
library(metastat) 
library(tidyverse)

#------------------loading data-------------------------------# 
load("HRB_sub_species_metadata.rds")
load("Jie_ASCVD_species_metadata.rds")

#------------------loading function -------------------------------# 
multiCompare <- function(dat, grp, type = "phenotype", covariates = NULL) {
  
  ## ---------- 0. preprocess ----------
  grp <- na.omit(grp)
  inte <- intersect(rownames(dat), rownames(grp))
  dat   <- dat[inte, , drop = FALSE]
  grp   <- grp[inte, , drop = FALSE]
  
  if (min(dat, na.rm = TRUE) >= 0 & max(dat, na.rm = TRUE) < 1 &&
      all(vapply(dat, is.numeric, TRUE))) {
    type <- "profile"
  }
  cat("The input data type is: ", type, "\n")
  
  grp[, 1] <- as.factor(grp[, 1])
  grp.level     <- levels(grp[, 1])
  grp.level.n   <- length(grp.level)
  if (grp.level.n != 2) stop("Number of Group levels must be 2")
  
  ## ---------- 1. logit.OR ----------
  out.cn <- paste(rep(c("median", "mean", "SD", "mean_rank",
                        "occ_rate", "n"), each = 2),
                  rep(grp.level, 6), sep = ".")
  out.cn <- c(c("pvalue", "FDR", "Enrichment", "Effectsize"), out.cn,
              "logit.pvalue", "logit.coeff", "logit.coeff.SE", "logit.coeff.pval", "logit.OR")
  
  out <- matrix(NA, ncol(dat), length(out.cn))
  out <- as.data.frame(out)
  colnames(out) <- out.cn
  rownames(out) <- colnames(dat)
  

  grpName <- colnames(grp)[1]
  if (is.null(covariates)) {
    rhs <- grpName
  } else if (is.data.frame(covariates)) {
    covariates <- covariates[inte, , drop = FALSE]
    rhs <- paste0(grpName, "~", paste(colnames(covariates), collapse = " + "))
  } else if (is.character(covariates)) {
    rhs <- paste0(grpName, "~", paste(colnames(covariates), collapse = " + "))
  } else {
    stop("covariates must be NULL, a data.frame or a character formula fragment")
  }
  
 
  a <- cbind(dat, grp)
  if (is.data.frame(covariates)) a <- cbind(a, covariates[inte,])
  
  for (i in rownames(out)) {
    if (all(is.na(dat[, i]))) next
   
    f <- as.formula(paste0(i, " ~ ", grpName))
    out[i, "pvalue"] <- wilcox.test(f, data = a)$p.value
    out[i, "Effectsize"] <- rstatix::wilcox_effsize(f, data = a)$effsize
    
    out[i, 5:6] <- tapply(dat[, i], grp[, 1], function(y) median(y, na.rm = TRUE))
    out[i, 7:8] <- tapply(dat[, i], grp[, 1], function(y) mean(y, na.rm = TRUE))
    out[i, 9:10] <- tapply(dat[, i], grp[, 1], function(y) sd(y, na.rm = TRUE))
    tmp  <- dat[, i]
    tmp2 <- tmp[!is.na(tmp)]
    out[i, 11:12] <- tapply(rank(tmp2), grp[, 1][!is.na(tmp)], function(y) mean(y))
    out[i, 13:14] <- tapply(dat[, i], grp[, 1],
                            function(y) ifelse(type == "phenotype",
                                               sum(!is.na(y))/length(y),
                                               sum(y > 0)/length(y)))
    out[i, 15:16] <- tapply(dat[, i], grp[, 1],
                            function(y) ifelse(type == "phenotype",
                                               sum(!is.na(y)),
                                               sum(y > 0)))
    
  
    a[,i] <- log(ifelse(a[,i]==0, min(a[,i][a[,i]!=0])/2, a[,i]))
    logit.form <- as.formula(paste0(rhs, " + ", i))
    
    fit <- try(glm(logit.form, data = a, family = binomial()), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      coef.sum <- summary(fit)$coefficients
    
      out[i, "logit.pvalue"]   <- coef.sum[i, "Pr(>|z|)"]
      out[i, "logit.coeff"]    <- coef.sum[i, "Estimate"]
      out[i, "logit.coeff.SE"] <- coef.sum[i, "Std. Error"]
      out[i, "logit.coeff.pval"] <- out[i, "logit.pvalue"]
      out[i, "logit.OR"]       <- exp(coef.sum[i, "Estimate"])  
    }
  }
  
  out$FDR <- p.adjust(out$pvalue, method = "BH")
  out$Enrichment <- ifelse(out[, 11] > out[, 12], grp.level[1], grp.level[2])
  out$Enrichment[out$pvalue > 0.05 | is.na(out$FDR)] <- "NONE"
  a <- ifelse(out[, 11] > out[, 12], -1, 1)
  out$Effectsize <- sign(a) * out$Effectsize
  
  out
}

#--------------------CKD--------------------------------------------# 

CKD_metadata$egfr_g <- factor(CKD_metadata$egfr_g, levels = c("NGT", "CKD"))

comres_egfr <- multiCompare(dat = t(heb_spe_egfr)/100, grp = CKD_metadata[,"egfr_g",drop=F], 
                           covariates = egfr[,c("Age", "Sex", "Smoking", "Alcohol",  "Drug_NSAID", 
                                              "Drug_Antibiotics" , "Drug_Diabetes" , "Drug_Hypertension", 
                                              "Drug_Dyslipidemia", "veg_freq")], type = "meta")
comres_egfr$logit.fdr <- p.adjust(comres_egfr$logit.coeff.pval, method = "BH")

#--------------------T2D--------------------------------------------# 
T2D_metadata$T2D_g <- factor(T2D_metadata$T2D_g, levels = c("NGT", "T2D"))

comres_t2d <- multiCompare(dat = t(heb_spe_t2d)/100, grp = T2D_metadata[,"T2D_g",drop=F], 
                            covariates = T2D_metadata[,c("Age", "Sex", "Smoking", "Alcohol",  "Drug_NSAID", 
                                                 "Drug_Antibiotics" , "Drug_Diabetes" , "Drug_Hypertension", 
                                                 "Drug_Dyslipidemia","veg_freq")], type = "meta")

comres_t2d$logit.fdr <- p.adjust(comres_t2d$logit.coeff.pval, method = "BH")


#------------------- Jie_ASCVD --------------------------------------------# 
jie_phe_f <- Jie_metadata[colnames(Jie_spe_pro), ]
jie_phe_f$group <- as.factor(ifelse(jie_phe_f$label1 == 0, "Control", "ASCVD"))
jie_phe_f$group <- factor(jie_phe_f$group, levels = c("Control", "ASCVD"))
jie_phe_f$age <- jie_phe_f$`Age.(year)`
jie_phe_f$Gender <- as.factor(jie_phe_f$Gender)
jie_phe_f$BMI <- jie_phe_f$`Body.Mass.Index.(BMI)`

comres_ascvd <- multiCompare(dat = t(Jie_spe_pro)/100, grp = jie_phe_f[,"group",drop=F], 
                           covariates = jie_phe_f[,c("age", "Gender")], type = "meta")
comres_ascvd$logit.fdr <- p.adjust(comres_ascvd$logit.coeff.pval, method = "BH")

#---------------- output ----------------------------------------------
openxlsx::write.xlsx(list(comres_t2d, comres_egfr, comres_ascvd), 
file = "T2D_CKD_ASCVD_validation_com_0726_n359.xlsx", rowNames = T)

