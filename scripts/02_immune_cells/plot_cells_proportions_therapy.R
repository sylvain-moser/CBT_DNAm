library(dplyr)

#prefix_statgen=ifelse(file.exists("/psycl/g/mpsstatgen/symo/dummy_file.txt"),"/psycl/g/mpsstatgen/symo/anxiety_methylation/","/Users/sylvain_moser/psycl_statgen/symo/anxiety_methylation/")
plot_dir=snakemake@params[["plot_dir"]]
#meth_df=readRDS(paste0(prefix_statgen,"data_smoking/09_immune_cells/cell_types_df.RData"))
meth_df=readRDS(snakemake@input[["cell_types_df"]])
meth_df =meth_df %>% filter(Timepoint %in% c("T0","T4","E","K1"))
meth_df$NID=factor(meth_df$NID) 

meth_df=meth_df %>% mutate(mdNLR=Gran/(CD8T+CD4T+NK+Bcell))

########################## plotting the results ###################
library(ggplot2)
time_df=meth_df %>% mutate(Time=case_when(
  Timepoint=="T0" ~ 0,
  Timepoint=="T4" ~ 1,
  Timepoint=="E" ~ 2,
  Timepoint=="K1" ~ 3
))
plot_Gran=ggplot(time_df,aes(x=Time,y=Gran,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "Granulocytes",y="Granulocytes cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/Granulocytes_therapy.png"),plot_Gran,width = 8.2,height = 7.5)

plot_CD8T=ggplot(time_df,aes(x=Time,y=CD8T,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "CD8T",y="CD8T cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/CD8T_cells_therapy.png"),plot_CD8T,width = 8.2,height = 7.5)

plot_CD4T=ggplot(time_df,aes(x=Time,y=CD4T,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "CD4T",y="CD4T cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=snakemake@output[["CD4T_plot"]],plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(time_df,aes(x=Time,y=mdNLR,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "mdNLR",y="mdNLR cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_therapy.png"),plot_mdNLR,width = 8.2,height = 7.5)

## color-coded for the remission and response

pheno=read.csv(snakemake@input[["pheno_file"]],sep=";")
pheno_df=pheno %>%
  mutate(response=case_when(HAMAend/HAMA <=0.5 ~ 1,HAMAend/HAMA > 0.5 ~ 0)) %>% 
  mutate(remission=case_when(HAMAend <=7 ~1,HAMAend >7 ~0)) %>% 
  dplyr::select(c(NID,response,remission))

remission_df=left_join(time_df,pheno_df) %>% filter (!is.na(remission))

# remission 

plot_Gran=ggplot(remission_df,aes(x=Time,y=Gran,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "Granulocytes",y="Granulocytes cells proportion",colour="remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/Granulocytes_therapy_remission.png"),plot_Gran,width = 8.2,height = 7.5)

plot_CD8T=ggplot(remission_df,aes(x=Time,y=CD8T,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "CD8T",y="CD8T cells proportion",colour="remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/CD8T_cells_therapy_remission.png"),plot_CD8T,width = 8.2,height = 7.5)

plot_CD4T=ggplot(remission_df,aes(x=Time,y=CD4T,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "CD4T",y="CD4T cells proportion",colour="remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/CD4T_cells_therapy_remission.png"),plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(remission_df,aes(x=Time,y=mdNLR,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "mdNLR",y="mdNLR cells proportion",colour="remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_therapy_remission.png"),plot_mdNLR,width = 8.2,height = 7.5)

# response

plot_Gran=ggplot(remission_df,aes(x=Time,y=Gran,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "Granulocytes",y="Granulocytes cells proportion",colour="response")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/Granulocytes_therapy_response.png"),plot_Gran,width = 8.2,height = 7.5)

plot_CD8T=ggplot(remission_df,aes(x=Time,y=CD8T,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "CD8T",y="CD8T cells proportion",colour="response")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/CD8T_cells_therapy_response.png"),plot_CD8T,width = 8.2,height = 7.5)

plot_CD4T=ggplot(remission_df,aes(x=Time,y=CD4T,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "CD4T",y="CD4T cells proportion",colour="response")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/CD4T_cells_therapy_response.png"),plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(remission_df,aes(x=Time,y=mdNLR,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "mdNLR",y="mdNLR cells proportion",colour="response")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_therapy_response.png"),plot_mdNLR,width = 8.2,height = 7.5)

##### centered on zero ##### 
time_df2=meth_df  %>% mutate(Time=case_when(
  Timepoint=="T0" ~ 0,
  Timepoint=="T4" ~ 1,
  Timepoint=="E" ~ 2,
  Timepoint=="K1" ~ 3
)) %>% filter (! NID %in% c("MPIPSYKL_014888","MPIPSYKL_014916","MPIPSYKL_014864","MPIPSYKL_013391","MPIPSYKL_012183","MPIPSYKL_0012908")) # individuals for which the 4 values are not available

#this will return an empty dataframe in the case of the example_random_data, thereofre:
if (dim(time_df2)[2]==0){
time_df2=meth_df %>% mutate(Time=case_when(
  Timepoint=="T0" ~ 0,
  Timepoint=="T4" ~ 1,
  Timepoint=="E" ~ 2,
  Timepoint=="K1" ~ 3
))
}

centered_time_df=time_df2 %>% group_by(NID) %>% group_map(~mutate(.,Gran=Gran  - as.numeric(.$Gran[.$Timepoint=="T0"])),.keep=T) %>% bind_rows()
centered_time_df=centered_time_df %>% group_by(NID) %>% group_map(~mutate(.,CD8T=CD8T  - as.numeric(.$CD8T[.$Timepoint=="T0"])),.keep=T) %>% bind_rows()
centered_time_df=centered_time_df %>% group_by(NID) %>% group_map(~mutate(.,CD4T=CD4T - as.numeric(.$CD4T[.$Timepoint=="T0"])),.keep=T) %>% bind_rows()


plot_Gran=ggplot(centered_time_df,aes(x=Time,y=Gran,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "Granulocytes",y="Granulocytes cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/Granulocytes_therapy_centered.png"),plot_Gran,width = 8.2,height = 7.5)

plot_CD8T=ggplot(centered_time_df,aes(x=Time,y=CD8T,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "CD8T",y="CD8T cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/CD8T_cells_therapy_centered.png"),plot_CD8T,width = 8.2,height = 7.5)

plot_CD4T=ggplot(centered_time_df,aes(x=Time,y=CD4T,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "CD4T",y="CD4T cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/CD4T_cells_therapy_centered.png"),plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(centered_time_df,aes(x=Time,y=mdNLR,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2,3),labels=c("T0","T4","E","K1"))+
  labs(title = "mdNLR",y="mdNLR cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_therapy_centered.png"),plot_mdNLR,width = 8.2,height = 7.5)
