t_test_meth_diff=function(cg_island,t0df,t1df){
  tryCatch({
    pval=t.test(t0df[,cg_island],t1df[,cg_island],paired=T)$p.value
    pval_df=data.frame("island"=cg_island,"pval_t_test"=pval)
    return(pval_df)
  },error=function(err){
    pval_df=data.frame("island"=cg_island,"pval_t_test"=NA,"qval_t_test"=NA)
    return(pval_df)
  })
}

did_not_return_exposure=function(cpgs_list,residuals=TRUE,methylation_data,residuals_data){
  
  library(assertthat)
  library(dplyr)
  library(pbmcapply)
  library(data.table)
  
  if (residuals==TRUE){
    meth_df=readRDS(residuals_data)
    
  }else{
    load(methylation_data)
  }
  
  t0NID=meth_df %>% filter (Timepoint=="b1") %>% pull (NID)
  t1NID=meth_df %>% filter (Timepoint=="p24_1") %>% pull (NID)
  
  t0df= meth_df %>% filter (Timepoint=="b1") %>% filter(NID %in% intersect(t0NID,t1NID)) %>% arrange(NID)
  t1df= meth_df %>% filter (Timepoint=="p24_1")%>% filter(NID %in% intersect(t0NID,t1NID)) %>% arrange(NID)
  
  assert_that(are_equal(t0df$NID,t1df$NID))
  
  t_test_list=lapply(cpgs_list,t_test_meth_diff,t0df=t0df,t1df=t1df)
  
  t_test_df=as.data.frame(rbindlist(t_test_list))
  t_test_df$qval_t_test=p.adjust(t_test_df$pval_t_test,method = "fdr")
  
  return(t_test_df)
}


