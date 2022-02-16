library(dplyr)
library(pbmcapply)

load(snakemake@input[["processed_diff_meth_df"]])
final_df=final_processed_results
load(snakemake@input[["meth_df"]])

ranked_df=final_df %>% mutate(stat_rank=rank(best_lmm_pval)) %>% arrange(stat_rank)
ranked_df=ranked_df # not rank the whole df it would take too long. 

compute_biological_rank_therapy=function(cpg_island){
  means=meth_df %>% select(c(cpg_island,"Timepoint")) %>% group_by(Timepoint) %>% summarise(mean=mean(eval(parse(text=(cpg_island)))))
  beta0=means %>% filter(Timepoint=="T0") %>% pull(mean)
  beta1=means %>% filter(Timepoint=="K1") %>% pull(mean)
  delta_beta=beta1-beta0
  return(delta_beta)
}


compute_biological_rank_exposure=function(cpg_island){
  means=meth_df %>% select(c(cpg_island,"Timepoint")) %>% group_by(Timepoint) %>% summarise(mean=mean(eval(parse(text=(cpg_island)))))
  beta0_1=means %>% filter(Timepoint=="b1") %>% pull(mean)
  beta1_1=means %>% filter(Timepoint=="pe1") %>% pull(mean)
  delta_beta_1=beta1_1-beta0_1
  beta0_2=means %>% filter(Timepoint=="p24_1") %>% pull(mean)
  beta1_2=means %>% filter(Timepoint=="pe1") %>% pull(mean)
  delta_beta_2=beta1_2-beta0_2  
  return(max(abs(delta_beta_1),abs(delta_beta_2)))
}

### exposure need another function with difference between t0 and t1 or t1-t2
if (length(grep("therapy",snakemake@input[["processed_diff_meth_df"]]))==1){
  deltas=mclapply(ranked_df$island,compute_biological_rank_therapy,mc.cores = as.numeric(snakemake@params[["cores"]]))
} else if (length(grep("exposure",snakemake@input[["processed_diff_meth_df"]]))==1){
  deltas=mclapply(ranked_df$island,compute_biological_rank_exposure,mc.cores = as.numeric(snakemake@params[["cores"]]))
}


delta_df=data.frame("island"=ranked_df$island,"delta_beta"=unlist(deltas))
ranked_df=left_join(ranked_df,delta_df)
ranked_df=ranked_df %>% mutate (bio_rank=rank(-abs(delta_beta))) %>% mutate(sum_rank=stat_rank+bio_rank) # ranking -x ranks it in descending order. 
save(ranked_df,file=snakemake@output[["ranked_cpgs"]])

library(GenomicRanges)
library(methyAnalysis)
library(Repitools)
library(dplyr)

bed_core=readRDS(snakemake@input[["illumina_annot_file"]])
bed=right_join(bed_core,ranked_df %>% filter(island !="NA")) %>% dplyr::rename(chrom=chr) %>% dplyr::rename(end=stop) %>% dplyr::select (c("chrom","start","end","island","qval","delta_beta","stat_rank","bio_rank","sum_rank"))

colnames(bed)[1:3]=c("chrom","chromStart","chromEnd")
dmrs=makeGRangesFromDataFrame(bed,keep.extra.columns = T)

annotations=annotateDMRInfo(dmrs,"TxDb.Hsapiens.UCSC.hg19.knownGene",as.GRanges = FALSE,promoterRange = 2000)
annotations_df=annoGR2DF(annotations$sigDMRInfo)
save(annotations_df,file=snakemake@output[["ranked_annotated_cpgs"]])

write.csv(annotations_df,file = snakemake@output[["ranked_annotated_cpgs_csv"]],quote = F,row.names = F)
