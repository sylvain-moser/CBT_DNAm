library(dplyr)
library(haven)

# use the spss file from AE as a key to match sentrix ID to NID and timepoints

spss_file=read_sav(snakemake@input[["NID_timepoints_keys_file"]])
spss_file=as.data.frame(spss_file)
keys=spss_file %>% dplyr::select(NID,Sample_ID,timepoint,SentrixBarcode_A,SentrixPosition_A) %>% mutate(Timepoint=timepoint)


keys=keys %>% mutate(Sentrix_ID=paste0(SentrixBarcode_A,"_",SentrixPosition_A)) %>% filter(NID !="") %>% filter(Timepoint !="")

# add the cell counts
cellcount=read.table(snakemake@input[["cell_count_file"]])
NIDs=rownames(cellcount)
cellcount$Sentrix_ID=NIDs

cell_type_df=left_join(cellcount,keys %>% dplyr::select("NID","Sentrix_ID","Timepoint"))

# add sex and age
pheno_file=read.csv(snakemake@input[["pheno_file"]],sep=";")
cell_type_df=left_join(cell_type_df,pheno_file %>% dplyr::select(c("NID","age","Sex")))

cell_type_df=cell_type_df %>% dplyr::select(c("Sentrix_ID","NID","Timepoint","CD8T","CD4T","NK","Bcell","Mono","Gran","age","Sex"))

base::saveRDS(cell_type_df,file = snakemake@output[["cell_types_df"]])
