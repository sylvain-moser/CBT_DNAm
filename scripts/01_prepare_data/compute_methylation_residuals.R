library(assertthat)
library(dplyr)
library(pbmcapply)
library(data.table)


load(snakemake@input[["methylation_df"]])

compute_residuals=function(cg_island,df){
  print (cg_island)
  residuals=residuals(lm(as.formula(paste0(cg_island,"~ CD8T +CD4T + NK+ Bcell + Mono+ Gran")),data=df))
  return(residuals)
}

print("starting to compute the linear models")

residuals_list=mclapply(colnames(meth_df)[which(grepl("cg",colnames(meth_df)))],compute_residuals,df=meth_df,mc.cores=20)

residuals_df=as.data.frame(do.call("cbind",residuals_list))
colnames(residuals_df)=colnames(meth_df)[which(grepl("cg",colnames(meth_df)))]
residuals_df=cbind(residuals_df,meth_df$Timepoint)
colnames(residuals_df)[dim(residuals_df)[2]]="Timepoint"
residuals_df=cbind(residuals_df,meth_df$NID)
colnames(residuals_df)[dim(residuals_df)[2]]="NID"

saveRDS(residuals_df,file=snakemake@output[["methylation_residuals_df"]])