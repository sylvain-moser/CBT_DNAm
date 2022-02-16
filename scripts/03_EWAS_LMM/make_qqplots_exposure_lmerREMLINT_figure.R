library(ggplot2)
library(dplyr)
library(data.table)

full_df_list=list()
files=snakemake@input[["LMM_results"]]
for (idx in seq(1,length(files))){
  load(files[idx])
  full_df_list[[idx]]=results_df
}

full_df=rbindlist(full_df_list)
full_df=as.data.frame(full_df)

full_df=full_df %>% mutate_at(c("pval_lmm1","pval_lmm2","coef_time","coef_time2","AIC1","AIC2","pval_random_slope","coef_random_slope.Time","AIC_random_slope1","AIC_random_slope2","pval_random_slope2nd"),~ as.numeric(as.character(.x)))



make_qqplots_list=function(model,no_singular_only=FALSE,ymax=NULL){
  pval_col=paste0("pval_",model)
  if (no_singular_only==TRUE) {
    if (model=="lmm1" | model=="lmm2"){
      singular_col="singular"
    } else if (model=="random_slope"){
      singular_col="singular_random"
    } else if (model=="random_slope2nd"){
      singular_col="singular_random2nd"
    }
    full_df=full_df %>% filter (!!as.symbol(singular_col) ==FALSE)
  }
  uni_d=runif(n = length(full_df[,pval_col]),min=0,max=1)
  uni_d=uni_d[order(uni_d,decreasing = T)]
  uni_d=-log(uni_d,10)
  pvals=full_df[,pval_col][order(full_df[,pval_col],decreasing = T)]
  lamdaGC=sprintf("Lambda GC=%s",lambdaGC=round(median(qchisq(1-full_df[,pval_col],1),na.rm = T)/qchisq(0.5,1),4))
  if (is.null(ymax)){
    ymax=max(-log(pvals,10))
  }
  p=ggplot(data = NULL,aes(x=uni_d,y=-log(pvals,10))) + 
    geom_point() + 
    ylim(0,ymax) + 
    geom_line(aes(x=uni_d,y=uni_d),colour="black")+
    labs(x="-log10 Expected pvals",y="-log10 Observed pvals")+
    theme(text = element_text(size = 13))+
    annotate("text",label=lamdaGC,color="black",x=3,y=0,size=4)
  return(p)
}

library(gridExtra)

qqplot_list=lapply(c("lmm1","lmm2","random_slope","random_slope2nd"),make_qqplots_list,no_singular_only=FALSE)
n <- length(qqplot_list)
nCol <- 4
plot=grid.arrange(grobs=qqplot_list,ncol=nCol)
ggsave(plot,file=snakemake@output[["QQplots"]],width=10, height=3.5)

