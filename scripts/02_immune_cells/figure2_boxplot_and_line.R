library(dplyr)

# A : exposure 

meth_df=readRDS(snakemake@input[["cell_types_df"]])
meth_df =meth_df %>% filter(Timepoint %in% c("b1","pe1","p24_1"))
meth_df$NID=factor(meth_df$NID) 

meth_df=meth_df %>% mutate(mdNLR=Gran/(CD8T+CD4T+NK+Bcell))


library(ggplot2)
time_df=meth_df %>% mutate(Time=case_when(
  Timepoint=="b1" ~ 0,
  Timepoint=="pe1" ~ 1,
  Timepoint=="p24_1" ~ 4
))

library(tidyr)


pheno=read.csv(snakemake@input[["pheno_file"]],sep=";")
pheno_df=pheno %>%
  mutate(response=case_when(HAMAend/HAMA <=0.5 ~ "Yes",HAMAend/HAMA > 0.5 ~ "No")) %>% 
  mutate(remission=case_when(HAMAend <=7 ~"Yes",HAMAend >7 ~"No")) %>% 
  dplyr::select(c(NID,response,remission)) 



remission_df=left_join(time_df,pheno_df) %>% filter (!is.na(remission))

facet_df=remission_df %>% select(c("NID","Time","CD4T","NK","Gran","remission","mdNLR")) %>% rename("CD4 T-cells"="CD4T","Granulocytes"="Gran","NK Cells"="NK","GLR"="mdNLR") %>%  pivot_longer(c("CD4 T-cells","NK Cells","Granulocytes","GLR"),names_to = "cell_type",values_to = "proportion")
# make a dataframe with the p-values:

text=data.frame("cell_type"=c("CD4 T-cells","NK Cells","Granulocytes","GLR"),"y"=c(0.4,0.2,0.8,4),"label"=c("adj.p=7.4e-4","adj.p=2.8e-5","adj.p=4.3e-5","adj.p=7.4e-4"),"Time"=c(1,1,1,1))

exposure_plot=ggplot(facet_df,aes(x=Time,y=proportion,group=Time))+
  geom_boxplot(outlier.size=1,width=1)+
  geom_line(aes(x=Time,y=proportion,group=NID,color=as.factor(remission)),linetype="11",size=0.5,alpha=0.5)+
  geom_point(aes(x=Time,y=proportion,color=as.factor(remission)),size=1,alpha=0.5)+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,1,2.5,4),labels=c("0","1","//","24"))+
  labs(y="Cell-type proportion",colour="Remission")+
  theme(legend.position = "right",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 10),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  geom_text(data = text,aes(x=1.5, y=y, label=label),
            color="black",
            size=3)+
  facet_wrap(~cell_type,ncol=4,scales = "free")

# B : therapy

meth_df=readRDS(snakemake@input[["cell_types_df"]])
meth_df =meth_df %>% filter(Timepoint %in% c("T0","E"))
meth_df$NID=factor(meth_df$NID) 

meth_df=meth_df %>% mutate(mdNLR=Gran/(CD8T+CD4T+NK+Bcell))

pheno=read.csv(snakemake@input[["pheno_file"]],sep=";")
pheno_df=pheno %>%
  mutate(response=case_when(HAMAend/HAMA <=0.5 ~ 1,HAMAend/HAMA > 0.5 ~ 0)) %>% 
  mutate(remission=case_when(HAMAend <=7 ~1,HAMAend >7 ~0)) %>% 
  select(c(NID,response,remission,Sex,age))

# compute LMM for statistical difference testing 

library(lmerTest)
  lmm_list=lapply(c("CD8T","CD4T","NK","Bcell","Mono","Gran","mdNLR"),FUN=function(cell_type){
  time_df=meth_df %>% select(NID,Timepoint,cell_type,age,Sex)
  remission_df=left_join(time_df,pheno_df)
  lmm=lmer(as.formula(paste0(cell_type,"~ Timepoint * remission + age + Sex + (1| NID)")),data=remission_df)
  anova_res=anova(lmm,type = 3)
  cell_df=data.frame("cell_type"=cell_type,"pvalue_time"=anova_res["Timepoint","Pr(>F)"],"pvalue_interaction"=anova_res["Timepoint:remission","Pr(>F)"],"pvalue_remission"=anova_res["remission","Pr(>F)"])
  return(cell_df)
})

lmm_df=data.frame(data.table::rbindlist(lmm_list))
lmm_df$qvalue_time=p.adjust(lmm_df$pvalue_time,method = "fdr")
lmm_df$qvalue_interaction=p.adjust(lmm_df$pvalue_interaction,method = "fdr")
lmm_df$qvalue_remission=p.adjust(lmm_df$pvalue_remission,method = "fdr")

# prepare the stat df for plotting the pvalues: 
stat_list=list()
for (cell_type in c("CD4T","CD8T","Bcell","Gran","mdNLR")){
  cell_df=data.frame("cell_type"=cell_type,"group1"="T0","group2"="E","pval"=paste0("adj.p=",format(lmm_df[which(lmm_df$cell_type==cell_type),"qvalue_time"],scientific=T,digits = 2)))
  stat_list[[cell_type]]=cell_df
}
stat_df=data.table::rbindlist(stat_list)
stat_df=stat_df%>% 
  mutate (y.position=case_when(
    cell_type=="CD4T" ~0.3,
    cell_type=="CD8T" ~ 0.175,
    cell_type=="Gran" ~ 0.825,
    cell_type=="Bcell" ~ 0.115,
    cell_type=="mdNLR" ~ 6.5
  )) %>%  mutate(cell_type=case_when(
  cell_type=="CD4T" ~"CD4 T-cells",
  cell_type=="CD8T" ~ "CD8 T-cells",
  cell_type=="Gran" ~ "Granulocytes",
  cell_type=="Bcell" ~ "B-cells",
  cell_type=="mdNLR" ~ "GLR")) 
stat_df$Timepoint=c(1,1,1,1,1)
# plotting 
library(ggpubr)
library(ggplot2)

remission_df=left_join(meth_df,pheno_df) %>% filter (!is.na(remission))
facet_df=remission_df %>% select(c("NID","Timepoint","CD4T","CD8T","Gran","remission","mdNLR","Bcell")) %>% rename("CD4 T-cells"="CD4T","Granulocytes"="Gran","CD8 T-cells"="CD8T","GLR"="mdNLR","B-cells"="Bcell") %>%  pivot_longer(c("CD4 T-cells","B-cells","Granulocytes","GLR","CD8 T-cells"),names_to = "cell_type",values_to = "proportion")
facet_df$Timepoint=factor(facet_df$Timepoint,levels=c("T0","E"),ordered = T)

therapy_plot=ggplot(facet_df,aes(x=Timepoint,y=proportion,group=Timepoint))+
  geom_boxplot(outlier.size=1,width=1)+
  geom_line(aes(x=Timepoint,y=proportion,group=NID),linetype="11",size=0.5,alpha=0.5)+
  geom_point(aes(x=Timepoint,y=proportion),size=1,alpha=0.5)+
  #scale_x_continuous(name="Hours after exposure",breaks=c(0,1),labels=c("0","1"))+
  labs(y="Cell-type proportion")+
  theme(legend.position = "right",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 10),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  geom_text(data = stat_df,aes(x=1.5, y=y.position, label=pval),
            color="black",
            size=3)+
  facet_wrap(~cell_type,ncol=5,scales = "free")
  
# put the figure together: 

figure2=ggarrange(exposure_plot,therapy_plot,labels=c("A","B"),nrow=2,common.legend = F)
ggsave(snakemake@output[["figure2_path"]],figure2,width = 180,height = 120,units = "mm",dpi = 600)




