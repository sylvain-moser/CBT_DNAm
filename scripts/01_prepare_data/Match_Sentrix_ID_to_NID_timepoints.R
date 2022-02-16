library(dplyr)
library(haven)

base::load(snakemake@input[["methylation_file"]])
df=data.frame(t(beta_value))

# use the spss file from AE as a key to match sentrix ID to NID and timepoints

spss_file=read_sav(snakemake@input[["NID_timepoints_keys_file"]])
spss_file=as.data.frame(spss_file)
keys=spss_file %>% dplyr::select(NID,Sample_ID,timepoint,SentrixBarcode_A,SentrixPosition_A) %>% mutate(Timepoint=timepoint)


keys=keys %>% mutate(Sentrix_ID=paste0(SentrixBarcode_A,"_",SentrixPosition_A)) %>% filter(NID !="") %>% filter(Timepoint !="")


# add the NID and Timepoints to the methylation data 
df$Sentrix_ID=rownames(df)
df = df %>% dplyr::select(Sentrix_ID,everything())
meth_df=left_join(df,keys %>% dplyr::select("NID","Sentrix_ID","Timepoint"))
# add the cell counts
cellcount=read.table(snakemake@input[["cell_count_file"]])
NIDs=rownames(cellcount)
cellcount=cellcount %>% mutate(across(everything(),~scale(.x))) #scaling the counts is better for lme4
cellcount$Sentrix_ID=NIDs

meth_df=left_join(meth_df,cellcount)

# add sex and age and smoking:
pheno_file=read.csv(snakemake@input[["pheno_file"]],sep=";")
meth_df=left_join(meth_df,pheno_file %>% dplyr::select(c("NID","age","Sex","Zigis")))

base::save(meth_df,file = snakemake@output[["methylation_df"]])

library(feather)

write_feather(meth_df,path = snakemake@output[["methylation_df_feather"]])
