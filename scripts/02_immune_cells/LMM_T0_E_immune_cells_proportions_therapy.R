library(data.table)
library(dplyr)

cell_types_df=readRDS(snakemake@input[["cell_types_df"]])

cell_types_df =cell_types_df %>% filter(Timepoint %in% c("T0","E"))
cell_types_df$NID=factor(cell_types_df$NID) 

cell_types_df=cell_types_df %>% mutate(mdNLR=Gran/(CD8T+CD4T+NK+Bcell))

pheno=read.csv(snakemake@input[["pheno_file"]],sep=";")
pheno_df=pheno %>%
  mutate(response=case_when(HAMAend/HAMA <=0.5 ~ 1,HAMAend/HAMA > 0.5 ~ 0)) %>% 
  mutate(remission=case_when(HAMAend <=7 ~1,HAMAend >7 ~0)) %>% 
  select(c(NID,response,remission,Sex,age))



library(ggplot2)

plot_list=lapply(c("CD8T","CD4T","NK","Bcell","Mono","Gran"),FUN=function(cell_type){
  time_df=cell_types_df %>% select(NID,Timepoint,cell_type)
  remission_df=left_join(time_df,pheno_df)
  remission_df$Timepoint=factor(remission_df$Timepoint,levels=c("T0","E"),ordered = T)
  plot_cell_type=ggplot(data=remission_df,aes(x=Timepoint,y=!!sym(cell_type),fill=Timepoint))+
    scale_fill_manual(values=c("#00BFC4","#F8766D"))+
    theme(legend.position = "none",panel.background = element_blank(),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
    geom_boxplot(outlier.size=1)+
    geom_line(aes(x=Timepoint,y=!!sym(cell_type),group=NID),colour="grey30", linetype="11",size=0.5)+
    geom_point(aes(x=Timepoint,y=!!sym(cell_type)),color="black",size=1)+
    facet_wrap(~ remission,labeller = label_both)
    return(plot_cell_type)
})

# for the mdNLR we need to zoom to exclude outliers 
cell_type="mdNLR"
time_df=cell_types_df %>% select(NID,Timepoint,cell_type)
remission_df=left_join(time_df,pheno_df)
remission_df$Timepoint=factor(remission_df$Timepoint,levels=c("T0","E"),ordered = T)
plot_cell_type=ggplot(data=remission_df,aes(x=Timepoint,y=!!sym(cell_type),fill=Timepoint))+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  theme(legend.position = "none",panel.background = element_blank(),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  geom_boxplot(outlier.size=1)+
  geom_line(aes(x=Timepoint,y=!!sym(cell_type),group=NID),colour="grey30", linetype="11",size=0.5)+
  geom_point(aes(x=Timepoint,y=!!sym(cell_type)),color="black",size=1)+
  facet_wrap(~ remission,labeller=label_both)
plot_list[[7]]=plot_cell_type

LMM_T0_E=gridExtra::grid.arrange(grobs=plot_list)

ggsave(filename = snakemake@output[["cell_types_results_plots"]],LMM_T0_E,dpi=600 ,width = 180,height=120 ,units = "mm")


library(lmerTest)
lmm_list=lapply(c("CD8T","CD4T","NK","Bcell","Mono","Gran","mdNLR"),FUN=function(cell_type){
  time_df=cell_types_df %>% select(NID,Timepoint,cell_type,age,Sex)
  remission_df=left_join(time_df,pheno_df)
  lmm=lmer(as.formula(paste0(cell_type,"~ Timepoint * remission + age + Sex + (1| NID)")),data=remission_df)
  anova_res=anova(lmm,type = 3)
  cell_df=data.frame("cell_type"=cell_type,"pvalue_time"=anova_res["Timepoint","Pr(>F)"],"pvalue_interaction"=anova_res["Timepoint:remission","Pr(>F)"],"pvalue_remission"=anova_res["remission","Pr(>F)"])
  return(cell_df)
})

lmm_df=rbindlist(lmm_list)
lmm_df$qvalue_time=p.adjust(lmm_df$pvalue_time,method = "fdr")
lmm_df$qvalue_interaction=p.adjust(lmm_df$pvalue_interaction,method = "fdr")
lmm_df$qvalue_remission=p.adjust(lmm_df$pvalue_remission,method = "fdr")


write.table(lmm_df,snakemake@output[["cell_types_results_table"]],row.names = F,quote = F)



