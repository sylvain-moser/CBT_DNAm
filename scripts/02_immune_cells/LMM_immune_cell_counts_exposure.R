compute_time_LMM=function(cell_type,dataset){
  f1_full=as.formula(paste0(cell_type,"~ Time + (1|NID)+ age + Sex"))
  lmm_1_full=lmer(f1_full,data=dataset,REML = TRUE)
  
  f2_full=as.formula(paste0(cell_type,"~ Time +I(Time^2) + (1|NID)+age + Sex"))
  lmm_2_full=lmer(f2_full,data=dataset,REML = TRUE)

  return (c("lmm_1_full"=lmm_1_full,"lmm_2_full"=lmm_2_full))
}

compute_time_random_slope_LMM=function(cell_type,dataset){
  f1_random_slope=as.formula(paste0(cell_type,"~ Time + (Time|NID)+ (1|NID)+ age + Sex"))
  lmm_1_random_slope=lmer(f1_random_slope,data=dataset,REML = TRUE)
  
  pval=anova(lmm_1_random_slope)["Time","Pr(>F)"]
  coef=fixef(lmm_1_random_slope)["Time"]
  AIC=extractAIC(lmm_1_random_slope)
  singular=isSingular(lmm_1_random_slope,tol=1e-5) # attention tolerance here is 1e-5 default whereas lmerControl has 1e-4 by default. They might not agree
  
  return(c("pval_random_slope"=pval,"coef_random_slope"=coef,"AIC_random_slope"=AIC,"singular_random"=as.character(singular),"convergence_random"=lmm_1_random_slope@optinfo$conv$opt)) #need to make character vector to report singular as T or F 
}

compute_time2nd_random_slope_LMM=function(cell_type,dataset){
  f2_random_slope=as.formula(paste0(cell_type,"~ I(Time^2) + (I(Time^2)|NID)+(1|NID)+age + Sex"))
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

compute_lmms_cell_counts=function(cell_type,cell_types_df){
  tryCatch({
  cell_type_df=cell_types_df %>% dplyr::select("NID","Timepoint","age","Sex",cell_type)
  cell_type_df=cell_type_df %>% mutate(Time=case_when(
    Timepoint=="b1" ~ 0,
    Timepoint=="pe1" ~ 1,
    Timepoint=="p24_1" ~ 2
  ))
  cell_type_df=cell_type_df %>% mutate(!!sym(cell_type) := as.numeric(!!sym(cell_type)))
  library(lmerTest)
  models=compute_time_LMM(cell_type,dataset = cell_type_df)
  diff_reg=compare_models(models)
  random_slope_results=compute_time_random_slope_LMM(cell_type,dataset = cell_type_df)
  random_slope2_results=compute_time2nd_random_slope_LMM(cell_type,dataset = cell_type_df)
  diff_reg=append(diff_reg,random_slope_results)
  diff_reg=append(diff_reg,random_slope2_results)
  diff_reg["cell_type"]=cell_type
  return(diff_reg)
  }, error=function(err) {
    diff_reg=c("best_model"=NA,"pval_lmm1"=NA,"pval_lmm2"=NA,"coef_time"=NA,"coef_time2"=NA,"AIC1"=NA,"AIC2"=NA,"singular"=NA,"convergence"=NA,"pval_random_slope"=NA,"coef_random_slope.Time"=NA,"AIC_random_slope1"=NA,"AIC_random_slope2"=NA,"singular_random"=NA,"convergence_random"=NA,"pval_random_slope2nd"=NA,"coef_random_slope2nd.I.Time.2."=NA,"AIC_random_slope2nd1"=NA,"AIC_random_slope2nd2"=NA,"singular_random2nd"=NA,"convergence_random2nd"=NA,"cell_type"=cell_type)
    return (diff_reg)
  })
}


####################### main ########################################################
library(data.table)
library(dplyr)

cell_types_df=readRDS(snakemake@input[["cell_types_df"]])
cell_types_df =cell_types_df %>% filter(Timepoint %in% c("b1","pe1","p24_1"))
cell_types_df$NID=factor(cell_types_df$NID) 

cell_types_df=cell_types_df %>% mutate(mdNLR=Gran/(CD8T+CD4T+NK+Bcell))


results=lapply(c("CD8T","CD4T","NK","Bcell","Mono","Gran","mdNLR"),compute_lmms_cell_counts,cell_types_df=cell_types_df)
list_df=lapply(results, as.data.frame.list)
results_df=data.table::rbindlist(list_df)
results_df=as.data.frame(results_df)

############################## pval correction 
full_df=results_df
# transform factor columns (introduced by data.table) into numeric:
full_df=full_df %>% mutate_at(c("pval_lmm1","pval_lmm2","coef_time","coef_time2","AIC1","AIC2","pval_random_slope","coef_random_slope.Time","AIC_random_slope1","AIC_random_slope2","pval_random_slope2nd","coef_random_slope2nd.I.Time.2.","AIC_random_slope2nd1","AIC_random_slope2nd2"),~ as.numeric(as.character(.x)))
full_df=full_df %>% mutate_at(c("best_model","cell_type"),~ as.character(.x))

# 1)  FDR correct the pvalues
library(tidyr)

pvals=full_df %>% select(c("cell_type","pval_lmm1","pval_lmm2","pval_random_slope","pval_random_slope2nd"))
pvals_long=pvals %>%
  pivot_longer(!cell_type,names_to = "model",values_to = "pval")
pvals_long$qval=p.adjust(pvals_long$pval,method = "fdr")
pvals_wide=pvals_long %>%
  select(c("cell_type","model","qval")) %>%
  pivot_wider(names_from = model,values_from = qval) %>%
  rename(qval_lmm1=pval_lmm1) %>% 
  rename(qval_lmm2=pval_lmm2) %>% 
  rename(qval_random_slope=pval_random_slope) %>% 
  rename(qval_random_slope2nd=pval_random_slope2nd)

full_df=left_join(full_df,pvals_wide)

check_model_ok=function(line,AIC_colname){
  if (AIC_colname=="AIC2"){
    return (isTRUE(line$convergence==0))
  } else if (AIC_colname=="AIC_random_slope2"){
    return (isTRUE(line$convergence_random==0))
  } else if (AIC_colname=="AIC_random_slope2nd2"){
    return (isTRUE(line$convergence_random2nd==0))
  }
}


get_best_model_name=function(AIC_colname){
  if (AIC_colname=="AIC2"){
    return ("random_intercept")
  } else if (AIC_colname=="AIC_random_slope2"){
    return ("random_slope")
  } else if (AIC_colname=="AIC_random_slope2nd2"){
    return ("random_slope2nd")
  }
} 


find_best_model=function(cell_type_name){
  cell_type_df=full_df %>% filter(cell_type==cell_type_name)
  AIC_order=names(sort(cell_type_df %>% select("AIC2","AIC_random_slope2","AIC_random_slope2nd2")))
  if (check_model_ok(cell_type_df,AIC_order[1])==TRUE){
    best_model=get_best_model_name(AIC_order[1])
  } else{
    if (check_model_ok(cell_type_df,AIC_order[2])==TRUE){
      best_model=get_best_model_name(AIC_order[2])
    } else { 
      if (check_model_ok(cell_type_df,AIC_order[3])==TRUE){
        best_model=get_best_model_name(AIC_order[3])
      }else{
        best_model=NA
      }
    }
  }
  return(best_model)
}


full_df=full_df %>% rename(best_intercept_model=best_model)    
full_df$best_model=sapply(full_df$cell_type,FUN = find_best_model)

process_results1=full_df %>% filter (best_model=="random_intercept") %>% 
  mutate(best_lmm_pval=case_when(best_intercept_model=="lmm_1" ~ pval_lmm1,best_intercept_model=="lmm_2" ~ pval_lmm2)) %>%
  mutate(qval=case_when(best_intercept_model=="lmm_1" ~ qval_lmm1,best_intercept_model=="lmm_2" ~ qval_lmm2)) %>%
  select (c("cell_type","best_model","coef_time","coef_time2","singular","convergence","best_lmm_pval","qval"))

process_results2=full_df %>% filter(best_model=="random_slope") %>% mutate(coef_time=coef_random_slope.Time) %>% mutate(coef_time2=0) %>% select(c("cell_type","best_model","coef_time","coef_time2","singular_random","convergence_random","pval_random_slope","qval_random_slope"))
process_results3=full_df %>% filter(best_model=="random_slope2nd") %>% mutate(coef_time=0) %>% mutate(coef_time2=coef_random_slope2nd.I.Time.2.) %>% select(c("cell_type","best_model","coef_time","coef_time2","singular_random2nd","convergence_random2nd","pval_random_slope2nd","qval_random_slope2nd"))


colnames(process_results2)=colnames(process_results1)
colnames(process_results3)=colnames(process_results1)

final_processed_results=rbind(process_results1,process_results2)
final_processed_results=rbind(final_processed_results,process_results3)
saveRDS(final_processed_results,file=snakemake@output[["cell_types_results_df"]])
write.table(final_processed_results,file=snakemake@output[["cell_types_results_table"]],row.names = F,quote = F)


## making the dataframe for the RF classification: 
immune_cells_wide=cell_types_df %>% select(NID,Timepoint,"CD4T","NK","Gran") %>%  pivot_wider(names_from = Timepoint,values_from = !c("NID","Timepoint"))
pheno=read.csv(snakemake@input[["pheno_file"]],sep=";")
pheno_df=pheno %>%
  mutate(response=case_when(HAMAend/HAMA <=0.5 ~ 1,HAMAend/HAMA > 0.5 ~ 0)) %>% 
  mutate(remission=case_when(HAMAend <=7 ~1,HAMAend >7 ~0)) %>% 
  select(c(NID,response,remission))
immune_cells_features_df=left_join(immune_cells_wide,pheno_df) %>% filter (!is.na(remission))
write.csv(immune_cells_features_df,snakemake@output[["immune_cells_features_df"]],row.names = F,quote = F)
