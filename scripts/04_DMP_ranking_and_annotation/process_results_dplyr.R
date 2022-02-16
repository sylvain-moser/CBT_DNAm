library(dplyr)
library(tidyr)
## get rid of the model which failed: but count them:
load(snakemake@input[["diff_meth_df"]])
full_df=full_df %>% filter (!is.na(best_model))
failed_models=full_df %>% filter (is.na(best_model)) %>% pull(island)

# transform factor columns (introduced by data.table) into numeric:
full_df=full_df %>% mutate_at(c("pval_lmm1","pval_lmm2","coef_time","coef_time2","AIC1","AIC2","pval_random_slope","coef_random_slope.Time","AIC_random_slope1","AIC_random_slope2","pval_random_slope2nd","coef_random_slope2nd.I.Time.2.","AIC_random_slope2nd1","AIC_random_slope2nd2"),~ as.numeric(as.character(.x)))
full_df=full_df %>% mutate_at(c("best_model","island"),~ as.character(.x))

# 1)  FDR correct the pvalues

pvals=full_df %>% select(c("island","pval_lmm1","pval_lmm2","pval_random_slope","pval_random_slope2nd"))
pvals_long=pvals %>%
  pivot_longer(!island,names_to = "model",values_to = "pval")
pvals_long$qval=p.adjust(pvals_long$pval,method = "fdr")
pvals_wide=pvals_long %>%
  select(c("island","model","qval")) %>%
  pivot_wider(names_from = model,values_from = qval) %>%
  rename(qval_lmm1=pval_lmm1) %>% 
  rename(qval_lmm2=pval_lmm2) %>% 
  rename(qval_random_slope=pval_random_slope) %>% 
  rename(qval_random_slope2nd=pval_random_slope2nd)

full_df=left_join(full_df,pvals_wide)

#2) compare AICs: select model with lower AIC:
# if this model converged: take that one and report 
# if this model failed to converge: check the other one
# if the second model doesnot fail to converge take the second
# if it does fail to converge take the first one (with lower AIC) -> write to fail dataframe


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


find_best_model=function(island_name){
  island_df=full_df %>% filter(island==island_name)
  AIC_order=names(sort(island_df %>% select("AIC2","AIC_random_slope2","AIC_random_slope2nd2")))
  if (check_model_ok(island_df,AIC_order[1])==TRUE){
    best_model=get_best_model_name(AIC_order[1])
  } else{
    if (check_model_ok(island_df,AIC_order[2])==TRUE){
      best_model=get_best_model_name(AIC_order[2])
    } else { 
      if (check_model_ok(island_df,AIC_order[3])==TRUE){
        best_model=get_best_model_name(AIC_order[3])
      }else{
        best_model=NA
      }
    }
  }
  return(best_model)
}


full_df=full_df %>% rename(best_intercept_model=best_model)    
full_df$best_model=sapply(full_df$island,FUN = find_best_model)

process_results1=full_df %>% filter (best_model=="random_intercept") %>% 
  mutate(best_lmm_pval=case_when(best_intercept_model=="lmm_1" ~ pval_lmm1,best_intercept_model=="lmm_2" ~ pval_lmm2)) %>%
  mutate(qval=case_when(best_intercept_model=="lmm_1" ~ qval_lmm1,best_intercept_model=="lmm_2" ~ qval_lmm2)) %>%
  select (c("island","best_model","coef_time","coef_time2","singular","convergence","best_lmm_pval","qval"))

process_results2=full_df %>% filter(best_model=="random_slope") %>% mutate(coef_time=coef_random_slope.Time) %>% mutate(coef_time2=0) %>% select(c("island","best_model","coef_time","coef_time2","singular_random","convergence_random","pval_random_slope","qval_random_slope"))
process_results3=full_df %>% filter(best_model=="random_slope2nd") %>% mutate(coef_time=0) %>% mutate(coef_time2=coef_random_slope2nd.I.Time.2.) %>% select(c("island","best_model","coef_time","coef_time2","singular_random2nd","convergence_random2nd","pval_random_slope2nd","qval_random_slope2nd"))

failed_models2=full_df %>% filter(is.na(best_model)) %>% pull (island)
failed_models=rbind(failed_models,failed_models2)
write.table(failed_models,file=snakemake@output[["failed_models"]])


colnames(process_results2)=colnames(process_results1)
colnames(process_results3)=colnames(process_results1)

final_processed_results=rbind(process_results1,process_results2)
final_processed_results=rbind(final_processed_results,process_results3)
save(final_processed_results,file=snakemake@output[["processed_diff_meth_df"]])

library(feather)

write_feather(final_processed_results,path = snakemake@output[["processed_diff_meth_df_feather"]])

                            