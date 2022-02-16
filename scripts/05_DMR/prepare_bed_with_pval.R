library(dplyr)

load(snakemake@input[["diff_meth_df"]])
pval_col=snakemake@params[["pval_column"]]

bed_core=readRDS(snakemake@input[["illumina_annot_file"]])
## use the random slope of the second order model as input for comb-p
bed=right_join(bed_core,final_processed_results %>% filter(island !="NA")) %>% dplyr::rename(chrom=chr) %>% dplyr::rename(end=stop) %>% select (c("chrom","start","end",pval_col)) %>%dplyr::rename(pval=pval_col)  
bed=bed %>% arrange_at(c("chrom","start"))
bed=bed %>% filter (!is.na(pval))
write.table(bed,file=snakemake@output[["regulated_island_bim_file"]],row.names = F,quote = F,sep="\t")



  