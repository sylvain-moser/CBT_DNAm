source(as.character(snakemake@params[["surrogate_variable_script"]]))

compute_time_LMM=function(cg_island,dataset,n_surogates){
  f1_full=as.formula(paste0(cg_island,"~ Time + (1|NID)+",paste0(c(sprintf("surogate_%s",seq(1,n_surogates))),collapse = "+"),"+CD8T +CD4T + NK+ Bcell + Mono+ Gran+ age + Sex + Zigis"))
  lmm_1_full=lmer(f1_full,data=dataset,REML = TRUE)
  
  f2_full=as.formula(paste0(cg_island,"~ Time +I(Time^2) + (1|NID)+",paste0(c(sprintf("surogate_%s",seq(1,n_surogates))),collapse = "+"),"+CD8T +CD4T + NK+ Bcell + Mono+ Gran+ age + Sex + Zigis"))
  lmm_2_full=lmer(f2_full,data=dataset,REML = TRUE)

  return (c("lmm_1_full"=lmm_1_full,"lmm_2_full"=lmm_2_full,"lmm_null"=lmm_null))
}

compute_time_random_slope_LMM=function(cg_island,dataset,n_surogates){
  f1_random_slope=as.formula(paste0(cg_island,"~ Time + (Time|NID)+ (1|NID)+",paste0(c(sprintf("surogate_%s",seq(1,n_surogates))),collapse = "+"),"+CD8T +CD4T + NK+ Bcell + Mono+ Gran+ age + Sex + Zigis"))
  lmm_1_random_slope=lmer(f1_random_slope,data=dataset,REML = TRUE)
  
  pval=anova(lmm_1_random_slope)["Time","Pr(>F)"]
  coef=fixef(lmm_1_random_slope)["Time"]
  AIC=extractAIC(lmm_1_random_slope)
  singular=isSingular(lmm_1_random_slope,tol=1e-5) # attention tolerance here is 1e-5 default whereas lmerControl has 1e-4 by default. They might not agree
  
  return(c("pval_random_slope"=pval,"coef_random_slope"=coef,"AIC_random_slope"=AIC,"singular_random"=as.character(singular),"convergence_random"=lmm_1_random_slope@optinfo$conv$opt)) #need to make character vector to report singular as T or F 
}

compute_time2nd_random_slope_LMM=function(cg_island,dataset,n_surogates){
  f2_random_slope=as.formula(paste0(cg_island,"~ I(Time^2) + (I(Time^2)|NID)+(1|NID)+",paste0(c(sprintf("surogate_%s",seq(1,n_surogates))),collapse = "+"),"+CD8T +CD4T + NK+ Bcell + Mono+ Gran + age + Sex + Zigis"))
  lmm_2_random_slope=lmer(f2_random_slope,data=dataset,REML = TRUE)
  
  pval=anova(lmm_2_random_slope)["I(Time^2)","Pr(>F)"]
  coef=fixef(lmm_2_random_slope)["I(Time^2)"]
  AIC=extractAIC(lmm_2_random_slope)
  singular=isSingular(lmm_2_random_slope,tol=1e-5) # attention tolerance here is 1e-5 default whereas lmerControl has 1e-4 by default. They might not agree
  
  return(c("pval_random_slope2nd"=pval,"coef_random_slope2nd"=coef,"AIC_random_slope2nd"=AIC,"singular_random2nd"=as.character(singular),"convergence_random2nd"=lmm_2_random_slope@optinfo$conv$opt)) #need to make character vector to report singular as T or F 
}

compare_models=function(models){
  LRT1_2=anova(models$lmm_1_full,models$lmm_2_full)$`Pr(>Chisq)`[2] 
  if (LRT1_2 > 0.05) { #model 1 is better
    chosen_model=models$lmm_1_full
    chosen_name="lmm_1"
    coefficient1=fixef(models$lmm_1_full)[["Time"]]
    coefficient2=0
    AIC=extractAIC(models$lmm_1_full)
    singular=isSingular(models$lmm_1_full)
  }
  else{ # model 2 is better
    chosen_model=models$lmm_2_full
    chosen_name="lmm_2" 
    coefficient1=fixef(models$lmm_2_full)[["Time"]]
    coefficient2=fixef(models$lmm_2_full)[["I(Time^2)"]]
    AIC=extractAIC(models$lmm_2_full)
    singular=isSingular(models$lmm_2_full)
  }
  return(c("best_model"=chosen_name,"pval_lmm1"=anova(models$lmm_1_full)["Time","Pr(>F)"],"pval_lmm2"=anova(models$lmm_2_full)["I(Time^2)","Pr(>F)"],"coef_time"=coefficient1,"coef_time2"=coefficient2,"AIC"=AIC,"singular"=singular,"convergence"=chosen_model@optinfo$conv$opt)) 
}

find_regulated_islands=function(island,sur_mat,surrogate_df_path=snakemake@params[["surrogate_df_path"]]){
  tryCatch({
  print(island)
  ## check if the surrogate variables for this island have been already computed:
  if(file.exists(paste0(snakemake@params[["surrogate_df_path"]],island,"_surrogates.RData"))==FALSE){
    print ("the surrogate variables were not found on disk and need to be computed")
    surogate_variables=compute_surrogates(dataset = meth_df,feature = island,surrogates_matrix =sur_mat,surrogate_df_path=surrogate_df_path)
  }
  else{
    print ("the surrogate variables were found on disk and do NOT need to be computed")
    load(paste0(snakemake@params[["surrogate_df_path"]],island,"_surrogates.RData"))
    surogate_variables=sv_df
  }
  #
  colnames(surogate_variables)=c(sprintf("surogate_%s",seq(1,dim(surogate_variables)[2])))
  island_df=meth_df %>% dplyr::select(c("NID","Timepoint",all_of(island),"CD8T","CD4T","NK","Bcell","Mono","Gran","age","Sex","Zigis"))
  island_df=cbind(island_df,surogate_variables)
  island_df=island_df %>% mutate(Time=case_when(
    Timepoint=="b1" ~ 0,
    Timepoint=="pe1" ~ 1,
    Timepoint=="p24_1" ~ 2
  ))
  island_df=island_df %>% mutate(!!sym(island) := rankTransPheno(!!sym(island),para_c = 3/8))
  library(lmerTest)
  models=compute_time_LMM(island,dataset = island_df,n_surogates = dim(surogate_variables)[2])
  diff_reg=compare_models(models)
  random_slope_results=compute_time_random_slope_LMM(island,dataset = island_df,n_surogates = dim(surogate_variables)[2])
  random_slope2_results=compute_time2nd_random_slope_LMM(island,dataset = island_df,n_surogates = dim(surogate_variables)[2])
  diff_reg=append(diff_reg,random_slope_results)
  diff_reg=append(diff_reg,random_slope2_results)
  diff_reg["island"]=island 
  return(diff_reg)
  }, error=function(err) {
    diff_reg=c("best_model"=NA,"pval_lmm1"=NA,"pval_lmm2"=NA,"coef_time"=NA,"coef_time2"=NA,"AIC1"=NA,"AIC2"=NA,"singular"=NA,"convergence"=NA,"pval_random_slope"=NA,"coef_random_slope.Time"=NA,"AIC_random_slope1"=NA,"AIC_random_slope2"=NA,"singular_random"=NA,"convergence_random"=NA,"pval_random_slope2nd"=NA,"coef_random_slope2nd.I.Time.2."=NA,"AIC_random_slope2nd1"=NA,"AIC_random_slope2nd2"=NA,"singular_random2nd"=NA,"convergence_random2nd"=NA,"island"=island)
    return (diff_reg)
  })
}


####################### main ########################################################
# for sbatch execution only
# args <- commandArgs(trailingOnly = TRUE)
# start_idx=as.numeric(args[1])

start_idx=as.numeric(snakemake@params[["start_idx"]])


library(dplyr)
library(data.table)
library(FRGEpistasis)

load(snakemake@input[["methylation_df"]])
meth_df =meth_df %>% filter(Timepoint %in% c("b1","pe1","p24_1"))
meth_df$NID=factor(meth_df$NID) # need to reset the levels if not there is more levels than actual obs -> crashes lme4
## check if the surrogate matrix files exists: 
if(file.exists(snakemake@params[["surrogate_variable_matrix"]])==FALSE){
  print ("the surrogate matrix was not fond on disc and need to be computed")
  surrogate_matrix=get_surrogates_matrix(meth_df[,2:425585],scale = F,n_var = 10000)
  save(surrogate_matrix,file=snakemake@params[["surrogate_variable_matrix"]])
} else{
load(snakemake@params[["surrogate_variable_matrix"]])
}
library(parallel)
library(pbmcapply)
if (start_idx==0){
  print("analyzing snps 2:25585")
  results=pbmcapply::pbmclapply(colnames(meth_df)[2:25585],find_regulated_islands,sur_mat=surrogate_matrix,mc.cores = 40)
  } else{
    print(sprintf("analyzing snps %s to %s",(start_idx-1)*25000+25586,start_idx*25000+25585))
    results=pbmcapply::pbmclapply(colnames(meth_df)[((start_idx-1)*25000+25586):(start_idx*25000+25585)],find_regulated_islands,sur_mat=surrogate_matrix,mc.cores = 40)
    }
list_df=lapply(results, as.data.frame.list)
results_df=data.table::rbindlist(list_df)
results_df=as.data.frame(results_df)
save(results_df,file=snakemake@output[["find_DMP_results"]])

# satisfy dummy condition:
write.table(data.frame("dummy"="dummy"),file = snakemake@output[["dummy_output"]])

