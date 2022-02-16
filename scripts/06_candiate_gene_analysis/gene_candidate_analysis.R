library(dplyr)
library(tidyr)

### define the genes of interest: 

candidates=read.table(snakemake@input[["candidates"]],h=T,stringsAsFactors = F,na.strings = NA,sep="\t")

cpgs_from_other_studies=unique(c(as.vector(candidates %>% filter (shimada!="") %>% pull (shimada)),
                          as.vector(candidates %>% filter (petersen!="") %>% pull (petersen)),
                          as.vector(candidates %>% filter (ziegler!="") %>% pull (ziegler))))

###
## get rid of the model which failed: but count them:
load(snakemake@input[["diff_meth_df"]])
full_df=full_df %>% filter (!is.na(best_model))
failed_models=full_df %>% filter (is.na(best_model)) %>% pull(island)

# transform factor columns (introduced by data.table) into numeric:
full_df=full_df %>% mutate_at(c("pval_lmm1","pval_lmm2","coef_time","coef_time2","AIC1","AIC2","pval_random_slope","coef_random_slope.Time","AIC_random_slope1","AIC_random_slope2","pval_random_slope2nd","coef_random_slope2nd.I.Time.2.","AIC_random_slope2nd1","AIC_random_slope2nd2"),~ as.numeric(as.character(.x)))
full_df=full_df %>% mutate_at(c("best_model","island"),~ as.character(.x))

# filter the dataset for CpGs annotated to gene of interest:

bed_core=readRDS(snakemake@input[["illumina_annot_file"]])
bed=right_join(bed_core,full_df %>% filter(island !="NA")) %>% dplyr::rename(chrom=chr) %>% dplyr::rename(end=stop) %>% dplyr::select (c("chrom","start","end","island"))

colnames(bed)[1:3]=c("chrom","chromStart","chromEnd")
dmrs=GenomicRanges::makeGRangesFromDataFrame(bed,keep.extra.columns = T)

annotations=methyAnalysis::annotateDMRInfo(dmrs,"TxDb.Hsapiens.UCSC.hg19.knownGene",as.GRanges = FALSE,promoterRange = 2000)
annotations_df=Repitools::annoGR2DF(annotations$sigDMRInfo)


### FDR correction and processing for cpgs_of_interest:

interest_df=full_df %>% filter (island %in% cpgs_from_other_studies)

pvals=interest_df %>% dplyr::select(c("island","pval_lmm1","pval_lmm2","pval_random_slope","pval_random_slope2nd"))
pvals_long=pvals %>%
  pivot_longer(!island,names_to = "model",values_to = "pval")
pvals_long$qval=p.adjust(pvals_long$pval,method = "fdr")
pvals_wide=pvals_long %>%
  dplyr::select(c("island","model","qval")) %>%
  pivot_wider(names_from = model,values_from = qval) %>%
  dplyr::rename(qval_lmm1=pval_lmm1) %>% 
  dplyr::rename(qval_lmm2=pval_lmm2) %>% 
  dplyr::rename(qval_random_slope=pval_random_slope) %>% 
  dplyr::rename(qval_random_slope2nd=pval_random_slope2nd)

interest_df=left_join(interest_df,pvals_wide)

#2) compare AICs: dplyr::select model with lower AIC:
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
  island_df=interest_df %>% filter(island==island_name)
  AIC_order=names(sort(island_df %>% dplyr::select("AIC2","AIC_random_slope2","AIC_random_slope2nd2")))
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


interest_df=interest_df %>% dplyr::rename(best_intercept_model=best_model)    
interest_df$best_model=sapply(interest_df$island,FUN = find_best_model)

process_results1=interest_df %>% filter (best_model=="random_intercept") %>% 
  mutate(best_lmm_pval=case_when(best_intercept_model=="lmm_1" ~ pval_lmm1,best_intercept_model=="lmm_2" ~ pval_lmm2)) %>%
  mutate(qval=case_when(best_intercept_model=="lmm_1" ~ qval_lmm1,best_intercept_model=="lmm_2" ~ qval_lmm2)) %>%
  dplyr::select (c("island","best_model","coef_time","coef_time2","singular","convergence","best_lmm_pval","qval"))

process_results2=interest_df %>% filter(best_model=="random_slope") %>% mutate(coef_time=coef_random_slope.Time) %>% mutate(coef_time2=0) %>% dplyr::select(c("island","best_model","coef_time","coef_time2","singular_random","convergence_random","pval_random_slope","qval_random_slope"))
process_results3=interest_df %>% filter(best_model=="random_slope2nd") %>% mutate(coef_time=0) %>% mutate(coef_time2=coef_random_slope2nd.I.Time.2.) %>% dplyr::select(c("island","best_model","coef_time","coef_time2","singular_random2nd","convergence_random2nd","pval_random_slope2nd","qval_random_slope2nd"))

failed_models2=interest_df %>% filter(is.na(best_model)) %>% pull (island)
failed_models=rbind(failed_models,failed_models2)


colnames(process_results2)=colnames(process_results1)
colnames(process_results3)=colnames(process_results1)

final_processed_results=rbind(process_results1,process_results2)
final_processed_results=rbind(final_processed_results,process_results3)

# annotate the CpGs back to the genes: 

annotated_final_processed_results=left_join(final_processed_results,annotations_df %>% dplyr::select(island,GeneSymbol)) %>% arrange(best_lmm_pval)
write.table(annotated_final_processed_results,snakemake@output[["candidate_analysis_results"]],row.names = F,quote = F)
