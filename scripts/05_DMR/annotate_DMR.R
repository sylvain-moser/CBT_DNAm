library(GenomicRanges)
library(methyAnalysis)
library(Repitools)
library(dplyr)

dmr_file=read.table(snakemake@input[["DMR_result_file"]],stringsAsFactors = F)
dmr_file_filt=dmr_file %>% filter (V5 >=as.numeric(snakemake@params[["n_probe_threshold"]]))
print(snakemake@params[["n_probe_threshold"]])
dmrs_bed=dmr_file_filt[,c(1:3,5,7)] # extract the sidak corrected pvalue in column 7 
colnames(dmrs_bed)=c("chrom","chromStart","chromEnd","n_probe","pval")
dmrs=makeGRangesFromDataFrame(dmrs_bed,keep.extra.columns=TRUE)

annotations=annotateDMRInfo(dmrs,"TxDb.Hsapiens.UCSC.hg19.knownGene",as.GRanges = FALSE,promoterRange = 2000)
annotations_df=annoGR2DF(annotations$sigDMRInfo) %>% filter (pval!=0)

genes=annotations_df %>% select (EntrezID)
genes=unique(genes)
annotated_DMRs=annotations_df #%>% filter (PROMOTER==TRUE)
write.table(annotated_DMRs,snakemake@output[["annotated_DMR"]],row.names = F,quote = F)
write.table(genes,snakemake@output[["DMR"]],row.names = F,col.names = F,quote = F)

genes=annotations_df %>% filter (PROMOTER==TRUE) %>% select (EntrezID)
annotated_DMRs=annotations_df %>% filter (PROMOTER==TRUE)
write.table(annotated_DMRs,snakemake@output[["annotated_DMR_near_promoter"]],row.names = F,quote = F)
write.table(genes,snakemake@output[["DMR_near_promoter"]],row.names = F,col.names = F,quote = F)
