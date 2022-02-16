prefix_statgen=ifelse(file.exists("/psycl/g/mpsstatgen/symo/dummy_file.txt"),"/psycl/g/mpsstatgen/symo/anxiety_methylation/","/Users/sylvain_moser/anxiety_methylation/")


library(assertthat)
library(dplyr)
library(pbmcapply)
library(data.table)

expression_df=readRDS(snakemake@input[["expression_df"]])

compute_residuals=function(exp_probe,df){
  print (exp_probe)
  residuals=residuals(lm(as.formula(paste0(exp_probe,"~ NEUT + TCEL + MONO + DEND + BCEL")),data=df))
  return(residuals)
}

print("starting to compute the linear models")

residuals_list=lapply(colnames(expression_df)[which(grepl("ILMN_",colnames(expression_df)))],compute_residuals,df=expression_df)

residuals_df=as.data.frame(do.call("cbind",residuals_list))
colnames(residuals_df)=colnames(expression_df)[which(grepl("ILMN_",colnames(expression_df)))]
residuals_df=cbind(residuals_df,expression_df$Timepoint)
colnames(residuals_df)[dim(residuals_df)[2]]="Timepoint"
residuals_df=cbind(residuals_df,expression_df$NID)
colnames(residuals_df)[dim(residuals_df)[2]]="NID"

saveRDS(residuals_df,file=snakemake@output[["expression_residuals_df"]])
