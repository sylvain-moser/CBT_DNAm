library(dplyr)

plot_dir=snakemake@params[["plot_dir"]]
meth_df=readRDS(snakemake@input[["cell_types_df"]])
meth_df =meth_df %>% filter(Timepoint %in% c("b1","pe1","p24_1"))
meth_df$NID=factor(meth_df$NID) 

meth_df=meth_df %>% mutate(mdNLR=Gran/(CD8T+CD4T+NK+Bcell))

########################## plotting the results ###################
library(ggplot2)
time_df=meth_df %>% mutate(Time=case_when(
  Timepoint=="b1" ~ 0,
  Timepoint=="pe1" ~ 1,
  Timepoint=="p24_1" ~ 2
))
plot_Gran=ggplot(time_df,aes(x=Time,y=Gran,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2),labels=c("0","1","2"))+
  labs(title = "Granulocytes",y="Granulocytes cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/Granulocytes_exposure.png"),plot_Gran,width = 8.2,height = 7.5)

plot_NK=ggplot(time_df,aes(x=Time,y=NK,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2),labels=c("0","1","2"))+
  labs(title = "Natural Killer cells",y="NK cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/NK_cells_exposure.png"),plot_NK,width = 8.2,height = 7.5)

plot_CD4T=ggplot(time_df,aes(x=Time,y=CD4T,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2),labels=c("0","1","2"))+
  labs(title = "CD4T",y="CD4T cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=snakemake@output[["CD4T_plot"]],plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(time_df,aes(x=Time,y=mdNLR,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2),labels=c("0","1","2"))+
  labs(title = "mdNLR",y="mdNLR cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_exposure.png"),plot_mdNLR,width = 8.2,height = 7.5)
################### plotting centered on zero #######



centered_time_df=time_df %>% group_by(NID) %>% group_map(~mutate(.,Gran=Gran  - as.numeric(.$Gran[.$Timepoint=="b1"])),.keep=T) %>% bind_rows()
centered_time_df=centered_time_df %>% group_by(NID) %>% group_map(~mutate(.,NK=NK  - as.numeric(.$NK[.$Timepoint=="b1"])),.keep=T) %>% bind_rows()
centered_time_df=centered_time_df %>% group_by(NID) %>% group_map(~mutate(.,CD4T=CD4T - as.numeric(.$CD4T[.$Timepoint=="b1"])),.keep=T) %>% bind_rows()

plot_Gran=ggplot(centered_time_df,aes(x=Time,y=Gran,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2),labels=c("0","1","2"))+
  labs(title = "Granulocytes",y="delta cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+ 
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/Granulocytes_exposure_centered.png"),width = 8.2,height = 7.5)

plot_NK=ggplot(centered_time_df,aes(x=Time,y=NK,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2),labels=c("0","1","2"))+
  labs(title = "Natural Killer cells",y="delta cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/NK_cells_exposure_centered.png"),width = 8.2,height = 7.5)

plot_CD4T=ggplot(centered_time_df,aes(x=Time,y=CD4T,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2),labels=c("0","1","2"))+
  labs(title = "CD4T",y="delta cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/CD4T_cells_exposure_centered.png"),plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(centered_time_df,aes(x=Time,y=mdNLR,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Timepoints",breaks=c(0,1,2),labels=c("0","1","2"))+
  labs(title = "mdNLR",y="delta cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.ticks=element_blank(),axis.line = element_line(colour ="darkgrey"))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_exposure_centered.png"),plot_mdNLR,width = 8.2,height = 7.5)

########################## plotting the results with the hours instead of timepoints  ###################
library(ggplot2)
time_df=meth_df %>% mutate(Time=case_when(
  Timepoint=="b1" ~ 0,
  Timepoint=="pe1" ~ 1,
  Timepoint=="p24_1" ~ 24
))
plot_Gran=ggplot(time_df,aes(x=Time,y=Gran,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Granulocytes",y="Granulocytes cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank())+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/Granulocytes_exposure_hours.png"),plot_Gran,width = 8.2,height = 7.5)

plot_NK=ggplot(time_df,aes(x=Time,y=NK,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Natural Killer cells",y="NK cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank())+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/NK_cells_exposure_hours.png"),plot_NK,width = 8.2,height = 7.5)

plot_CD4T=ggplot(time_df,aes(x=Time,y=CD4T,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "CD4T",y="CD4T cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank())+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/CD4T_cells_exposure_hours.png"),plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(time_df,aes(x=Time,y=mdNLR,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "mdNLR",y="mdNLR cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank())+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_exposure_hours.png"),plot_mdNLR,width = 8.2,height = 7.5)

## color-coded for the remission and response

pheno=read.csv(snakemake@input[["pheno_file"]],sep=";")
pheno_df=pheno %>%
  mutate(response=case_when(HAMAend/HAMA <=0.5 ~ 1,HAMAend/HAMA > 0.5 ~ 0)) %>% 
  mutate(remission=case_when(HAMAend <=7 ~1,HAMAend >7 ~0)) %>% 
  dplyr::select(c(NID,response,remission))

remission_df=left_join(time_df,pheno_df) %>% filter (!is.na(remission))

# remission 

plot_Gran=ggplot(remission_df,aes(x=Time,y=Gran,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Granulocytes",y="Granulocytes cells proportion",colour="remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/Granulocytes_exposure_hours_remission.png"),plot_Gran,width = 8.2,height = 7.5)

plot_NK=ggplot(remission_df,aes(x=Time,y=NK,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Natural Killer cells",y="NK cells proportion",colour="remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/NK_cells_exposure_hours_remission.png"),plot_NK,width = 8.2,height = 7.5)

plot_CD4T=ggplot(remission_df,aes(x=Time,y=CD4T,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="CD4T after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "CD4T",y="delta cells proportion",colour="remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/CD4T_cells_exposure_hours_remission.png"),plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(remission_df,aes(x=Time,y=mdNLR,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="mdNLR after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "mdNLR",y="delta cells proportion",colour="remission")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_exposure_hours_remission.png"),plot_mdNLR,width = 8.2,height = 7.5)

# response

plot_Gran=ggplot(remission_df,aes(x=Time,y=Gran,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Granulocytes",y="Granulocytes cells proportion",colour="response")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/Granulocytes_exposure_hours_response.png"),plot_Gran,width = 8.2,height = 7.5)

plot_NK=ggplot(remission_df,aes(x=Time,y=NK,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Natural Killer cells",y="NK cells proportion",colour="response")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/NK_cells_exposure_hours_response.png"),plot_NK,width = 8.2,height = 7.5)

plot_CD4T=ggplot(remission_df,aes(x=Time,y=CD4T,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "CD4T",y="CD4T cells proportion",colour="response")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/CD4T_cells_exposure_hours_response.png"),plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(remission_df,aes(x=Time,y=mdNLR,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "mdNLR",y="mdNLR cells proportion",colour="response")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_exposure_hours_response.png"),plot_mdNLR,width = 8.2,height = 7.5)



################### plotting centered on zero #######
centered_time_df=time_df %>% group_by(NID) %>% group_map(~mutate(.,Gran=Gran  - as.numeric(.$Gran[.$Timepoint=="b1"])),.keep=T) %>% bind_rows()
centered_time_df=centered_time_df %>% group_by(NID) %>% group_map(~mutate(.,NK=NK  - as.numeric(.$NK[.$Timepoint=="b1"])),.keep=T) %>% bind_rows()
centered_time_df=centered_time_df %>% group_by(NID) %>% group_map(~mutate(.,CD4T=CD4T - as.numeric(.$CD4T[.$Timepoint=="b1"])),.keep=T) %>% bind_rows()

plot_Gran=ggplot(centered_time_df,aes(x=Time,y=Gran,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Granulocytes",y="delta cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank())+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/Granulocytes_exposure_hours_centered.png"),plot_Gran,width = 8.2,height = 7.5)

plot_NK=ggplot(centered_time_df,aes(x=Time,y=NK,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Natural Killer cells",y="delta cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank())+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/NK_cells_exposure_hours_centered.png"),plot_NK,width = 8.2,height = 7.5)

plot_CD4T=ggplot(centered_time_df,aes(x=Time,y=CD4T,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "CD4T",y="delta cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank())+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/CD4T_cells_exposure_hours_centered.png"),plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(centered_time_df,aes(x=Time,y=mdNLR,group=NID,colour=NID))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "mdNLR",y="delta cells proportion")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank())+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=1))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_exposure_hours_centered.png"),plot_mdNLR,width = 8.2,height = 7.5)


################### color-coded for the remission and response status #######

remission_df=left_join(centered_time_df,pheno_df) %>% filter (!is.na(remission))

# remission 

plot_Gran=ggplot(remission_df,aes(x=Time,y=Gran,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Granulocytes",y="delta cells proportion")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/Granulocytes_exposure_hours_centered_remission.png"),plot_Gran,width = 8.2,height = 7.5)

plot_NK=ggplot(remission_df,aes(x=Time,y=NK,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Natural Killer cells",y="delta cells proportion")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/NK_cells_exposure_hours_centered_remission.png"),plot_NK,width = 8.2,height = 7.5)

plot_CD4T=ggplot(remission_df,aes(x=Time,y=CD4T,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "CD4T",y="delta cells proportion")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/CD4T_cells_exposure_hours_centered_remission.png"),plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(remission_df,aes(x=Time,y=mdNLR,group=NID,colour=as.factor(remission)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "mdNLR",y="delta cells proportion")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(remission)))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_exposure_hours_centered_remission.png"),plot_mdNLR,width = 8.2,height = 7.5)

# response

plot_Gran=ggplot(remission_df,aes(x=Time,y=Gran,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Granulocytes",y="delta cells proportion")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/Granulocytes_exposure_hours_centered_response.png"),plot_Gran,width = 8.2,height = 7.5)

plot_NK=ggplot(remission_df,aes(x=Time,y=NK,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "Natural Killer cells",y="delta cells proportion")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/NK_cells_exposure_hours_centered_response.png"),plot_NK,width = 8.2,height = 7.5)

plot_CD4T=ggplot(remission_df,aes(x=Time,y=CD4T,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "CD4T",y="delta cells proportion")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/CD4T_cells_exposure_hours_centered_response.png"),plot_CD4T,width = 8.2,height = 7.5)

plot_mdNLR=ggplot(remission_df,aes(x=Time,y=mdNLR,group=NID,colour=as.factor(response)))+geom_line()+
  scale_x_continuous(name="Hours after exposure",breaks=c(0,6,12,18,24),labels=c("0","6","12","18","24"))+
  labs(title = "mdNLR",y="delta cells proportion")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),text = element_text(size = 30),axis.line = element_line(colour ="darkgrey"),axis.ticks=element_blank(),legend.title=element_text(size = 20),legend.text = element_text(size = 20))+
  stat_summary(fun=mean,geom="line",lwd=2,aes(group=as.factor(response)))
ggsave(filename=paste0(plot_dir,"/mdNLR_cells_exposure_hours_centered_response.png"),plot_mdNLR,width = 8.2,height = 7.5)



