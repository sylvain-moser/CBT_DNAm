library(dplyr)
library(haven)

base::load(snakemake@input[["expression_file"]])
df=data.frame(t(cdata2a))

link1=read.csv(snakemake@input[["Sample_ID_keys"]],sep =  " ")
link2=read_sav(snakemake@input[["NID_timepoints_keys_file"]])

# match the Bead.Chip_Section to the NID used for methylation

df$Bead.Chip_Section=rownames(df)
df = df %>% dplyr::select(Bead.Chip_Section,everything())

expression_df=left_join(df,link1 %>% select("Bead.Chip_Section","Sample_ID"))

expression_df=right_join(link2 %>% select("Sample_ID","NID","timepoint"),expression_df) %>% mutate(Timepoint=timepoint) %>% select("NID","Timepoint","Bead.Chip_Section",everything())

# add the cell counts
cellcount=read.table(snakemake@input[["cell_composition_expression"]],h=T,stringsAsFactors = F)

cellcount=cellcount %>% mutate(across(-sample,~scale(.x))) #scaling the counts is better for lme4

expression_df=left_join(expression_df,cellcount,by=c("Bead.Chip_Section"="sample"))

# add sex and age and smoking
pheno_file=read.csv(snakemake@input[["pheno_file"]],sep=";")

expression_df=left_join(expression_df,pheno_file %>% dplyr::select(c("NID","age","Sex","Zigis")))

saveRDS(expression_df,snakemake@output[["expression_df"]])
