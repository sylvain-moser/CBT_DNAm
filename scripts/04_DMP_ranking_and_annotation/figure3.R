library(dplyr)
library(ggplot2)
library(ggplotify)
library(ggpubr)
library(gridExtra)

load(snakemake@input[["exposure_cpg_ranked_anotated"]])
exposure_annotated=annotations_df

load(snakemake@input[["exposure_cpg_ranked"]])
exposure_ranked=ranked_df

load(snakemake@input[["therapy_cpg_ranked_anotated"]])
therapy_annotated=annotations_df

load(snakemake@input[["therapy_cpg_ranked"]])
therapy_ranked=ranked_df


# plot 1 : model distribution
#exposure

exposure_annotated=left_join(exposure_ranked,exposure_annotated %>% dplyr::select("island","EntrezID","GeneSymbol"))
exposure_annotated=exposure_annotated %>% arrange(sum_rank)

best_exposure=exposure_annotated[1:100,]
random_intercept_1= best_exposure %>% filter (best_model=="random_intercept" & coef_time2==0)
random_intercept_2=best_exposure %>% filter (best_model=="random_intercept" & coef_time2!=0)
random_slope=best_exposure %>% filter (best_model=="random_slope")
random_slope2nd=best_exposure %>% filter (best_model=="random_slope2nd")

exposure_models=data.frame("Model"=c("Linear; Random Intercept","Non-linear; Random Intercept","Linear; Random slope","Non-linear; Random slope"),
                        "Selected"=c(length(random_intercept_1$island),length(random_intercept_2$island),length(random_slope$island),length(random_slope2nd$island)),"phase"=rep("exposure",4))

#therapy:
therapy_annotated=left_join(therapy_ranked,therapy_annotated %>% dplyr::select("island","EntrezID","GeneSymbol"))
therapy_annotated=therapy_annotated %>% arrange(sum_rank)

best_therapy=therapy_annotated[1:100,]
random_intercept_1= best_therapy %>% filter (best_model=="random_intercept" & coef_time2==0)
random_intercept_2=best_therapy %>% filter (best_model=="random_intercept" & coef_time2!=0)
random_slope=best_therapy %>% filter (best_model=="random_slope")
random_slope2nd=best_therapy %>% filter (best_model=="random_slope2nd")

therapy_models=data.frame("Model"=c("Linear; Random Intercept","Non-linear; Random Intercept","Linear; Random slope","Non-linear; Random slope"),
                           "Selected"=c(length(random_intercept_1$island),length(random_intercept_2$island),length(random_slope$island),length(random_slope2nd$island)),"phase"=rep("therapy",4))

models_df=rbind(exposure_models,therapy_models)


model_plot=ggplot(data=models_df,aes(x=phase,y=Selected,fill=Model))+
geom_bar(stat="identity",position="dodge")+
# scale_fill_discrete(guide=guide_legend(override.aes = list(size = 5)))+
theme(legend.position = "top",legend.direction = "vertical",text = element_text(size = 8.5),axis.title.x = element_blank(),,panel.background = element_blank(),axis.line = element_line(colour ="darkgrey"),
      legend.title=element_text(size = 6),legend.text = element_text(size = 6),legend.key.height=unit(3, "mm"))
#guides(color = guide_legend(override.aes = list(size=10)))



# plot 2: HTR3A all people 6 time points: 

meth_residuals=readRDS(snakemake@input[["methylation_residuals_df"]])

exposure_ID=meth_residuals %>% filter(Timepoint %in% c("b1","pe1","p24_1")) %>% pull(NID)
therapy_ID=meth_residuals %>% filter(Timepoint %in% c("T0","T4","E","K1")) %>% pull(NID)

time_df=meth_residuals %>% select(NID,Timepoint,cg01586609) %>% filter (NID %in% intersect(exposure_ID,therapy_ID)) %>% 
  mutate(Time=case_when(
    Timepoint=="T0" ~ 0,
    Timepoint=="T4" ~ 1,
    Timepoint=="b1" ~ 2,
    Timepoint=="pe1" ~ 3,
    Timepoint=="p24_1" ~ 4,
    Timepoint=="E" ~ 5,
    Timepoint=="K1" ~ 6,
  ))


plot_HTR3A_CBT_residuals=ggplot(time_df,aes(x=Time,y=cg01586609,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3,4,5,6),labels=c("T0","T4","BE","P1h","P24h","E","K"))+
  labs(title = "",y="Residuals methylytion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 8.6),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=1.5,aes(group=1))



# plot 3: remission centered 

centered_time_df=time_df %>% filter(Timepoint %in% c("b1","pe1","p24_1","E","K1")) %>% group_by(NID) %>% group_map(~mutate(.,cg01586609=cg01586609  - as.numeric(.$cg01586609[.$Timepoint=="b1"])),.keep=T) %>% bind_rows()
pheno=read.csv(snakemake@input[["pheno_file"]],sep=";")
pheno_df=pheno %>%
  mutate(response=case_when(HAMAend/HAMA <=0.5 ~ 1,HAMAend/HAMA > 0.5 ~ 0)) %>% 
  mutate(remission=case_when(HAMAend <=7 ~"Yes",HAMAend >7 ~"No")) %>% 
  select(c(NID,response,remission))

remission_df=left_join(centered_time_df,pheno_df)

plot_legend=ggplot(remission_df,aes(x=Time,y=cg01586609,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(2,3,4,5,6),labels=c("BE","P1h","P24h","E","K"))+
  labs(title = "",y="Delta residuals methylation",colour="Remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 8.6),axis.ticks=element_blank(),legend.title=element_text(size = 6),legend.text = element_text(size = 6),axis.line = element_line(colour ="darkgrey"),legend.position = "top")+
  stat_summary(fun=mean,geom="line",lwd=1.5,aes(group=as.factor(remission)))

common_legend=get_legend(plot_legend)

plot_HTR3A_CBT_residuals_centered_postint_remission=ggplot(remission_df,aes(x=Time,y=cg01586609,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(2,3,4,5,6),labels=c("BE","P1h","P24h","E","K"))+
  labs(title = "",y="Delta residuals methylation",colour="Remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 8.6),axis.ticks=element_blank(),legend.title=element_text(size = 8.6),legend.text = element_text(size = 8.6),axis.line = element_line(colour ="darkgrey"),legend.position = "top",legend.direction = "vertical")+
  stat_summary(fun=mean,geom="line",lwd=1.5,aes(group=as.factor(remission)))

# plot 4: 

exp_df=readRDS(snakemake@input[["expression_residuals_df"]])

exposure_ID=exp_df %>% filter(Timepoint %in% c("b1","pe1","p24_1")) %>% pull(NID)

time_df=exp_df %>% dplyr::select(NID,Timepoint,ILMN_1662070)  %>% 
  mutate(Time=case_when(
    Timepoint=="b1" ~0,
    Timepoint=="pe1" ~ 1,
    Timepoint=="p24_1" ~ 2
  ))


remission_df=left_join(time_df,pheno_df) %>% filter (!is.na(remission))

plot_ILMN_1662070_postint_remission=ggplot(remission_df,aes(x=Time,y=ILMN_1662070,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2),labels=c("BE","P1h","P24h"))+
  labs(title = "",y="Expression residuals",colour="remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 8.6),axis.ticks=element_blank(),legend.title=element_text(size = 8.6),legend.text = element_text(size = 8.6),axis.line = element_line(colour ="darkgrey"),legend.position = "top",legend.direction = "vertical")+
  stat_summary(fun=mean,geom="line",lwd=1.5,aes(group=as.factor(remission)))


# plot 5:

time_df=exp_df %>% dplyr::select(NID,Timepoint,ILMN_2371079)  %>% 
  mutate(Time=case_when(
    Timepoint=="b1" ~0,
    Timepoint=="pe1" ~ 1,
    Timepoint=="p24_1" ~ 2
  ))


remission_df=left_join(time_df,pheno_df) %>% filter (!is.na(remission))

plot_ILMN_2371079_postint_remission=ggplot(remission_df,aes(x=Time,y=ILMN_2371079,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2),labels=c("BE","P1h","P24h"))+
  labs(title = "",y="Expression residuals",colour="remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 8.6),axis.ticks=element_blank(),legend.title=element_text(size = 8.6),legend.text = element_text(size = 8.6),axis.line = element_line(colour ="darkgrey"),legend.position = "top",legend.direction = "vertical")+
  stat_summary(fun=mean,geom="line",lwd=1.5,aes(group=as.factor(remission)))

# plot6 : 

averaged_df=exp_df %>% dplyr::select(NID,Timepoint,ILMN_1662070,ILMN_2371079)  %>% 
  mutate(Time=case_when(
    Timepoint=="b1" ~0,
    Timepoint=="pe1" ~ 1,
    Timepoint=="p24_1" ~ 2
  )) %>% 
  mutate(average_exp=(ILMN_1662070+ILMN_2371079)/2)

averaged_remission_df=left_join(averaged_df,pheno_df) %>% filter (!is.na(remission))

plot_averaged_postint_remission=ggplot(averaged_remission_df,aes(x=Time,y=average_exp,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2),labels=c("BE","P1h","P24h"))+
  labs(title = "",y="Expression residuals",colour="remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 8.6),axis.ticks=element_blank(),legend.title=element_text(size = 8.6),legend.text = element_text(size = 8.6),axis.line = element_line(colour ="darkgrey"),legend.position = "top",legend.direction = "vertical")+
  stat_summary(fun=mean,geom="line",lwd=1.5,aes(group=as.factor(remission)))

# try to make figure with that 
figure3=grid.arrange(model_plot,plot_HTR3A_CBT_residuals,plot_HTR3A_CBT_residuals_centered_postint_remission+theme(legend.position ="hidden"),
             plot_ILMN_1662070_postint_remission+theme(legend.position = "hidden"),plot_ILMN_2371079_postint_remission + theme(legend.position = "hidden"),
             plot_averaged_postint_remission+theme(legend.position = "hidden"),bottom=common_legend$grobs[[1]],ncol=3)
ggsave(figure3,filename =snakemake@output[["figure3_path"]],dpi=600 ,width = 180,height=120 ,units = "mm")


