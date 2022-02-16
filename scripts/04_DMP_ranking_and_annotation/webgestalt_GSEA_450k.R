if (!"WebGestaltR" %in% (.packages())){
  library(WebGestaltR)
}

if (!"GenomicRanges" %in% (.packages())){
  library(GenomicRanges)
}
if (!"methyAnalysis" %in% (.packages())){
  library(methyAnalysis)
}


if (!exists("all_genes")){  
  pos=readRDS("/psycl/g/mpsstatgen/symo/PD_exposure/data/prepared_data/illumina_450k_annotations.RDS")
  pos=pos %>% dplyr::select(c("chr","start","stop"))
  colnames(pos)=c("chr","chromStart","chromEnd")
  GR=makeGRangesFromDataFrame(pos,keep.extra.columns=TRUE)
  all_genes=annotateDMRInfo(GR,"TxDb.Hsapiens.UCSC.hg19.knownGene",as.GRanges = FALSE,promoterRange = 2000)
  all_genes=as.data.frame(all_genes$sigDMRInfo)
}
  
webgestalt_ORA=function(interest_genes_EntrezID,referenceGene_EntrezID=all_genes$EntrezID,FDR_threshold=0.1){
  print (sprintf("ORA with Webgestalt using FDR %s",FDR_threshold))
  ORA=WebGestaltR(enrichMethod = "ORA",
                  enrichDatabase = c("geneontology_Biological_Process","pathway_KEGG","pathway_Reactome"),
                  interestGene = interest_genes_EntrezID ,
                  interestGeneType = "entrezgene",
                  referenceGene = referenceGene_EntrezID ,
                  referenceGeneType = "entrezgene",
                  isOutput = F,minNum = 3,
                  fdrThr=FDR_threshold)
  if (exists("ORA") & !is.null(ORA)) {
    if (dim(ORA)[1] < 10){
    print (ORA[,c("geneSet","description","enrichmentRatio","pValue","FDR")])
    } else{
      print(sprintf("There are %s significant Gene sets. Displaying only the %s most significant ones",dim(ORA)[1],10))
      print (ORA[1:10,c("geneSet","description","enrichmentRatio","pValue","FDR")])
    }
    #rm(ORA)
  } else{
      print ("No significant results...")
    }
}


webgestalt_GSEA=function(interest_genes_EntrezID,referenceGene_EntrezID=all_genes$EntrezID,FDR_threshold=0.1){
  print (sprintf("GSEA with Webgestalt using FDR %s",FDR_threshold))
  GSEA=WebGestaltR(enrichMethod = "GSEA",
                  enrichDatabase = c("geneontology_Biological_Process","pathway_KEGG","pathway_Reactome"),
                  interestGene = interest_genes_EntrezID ,
                  interestGeneType = "entrezgene",
                  referenceGene = referenceGene_EntrezID ,
                  referenceGeneType = "entrezgene",
                  isOutput = F,minNum = 3,
                  fdrThr=FDR_threshold)
  if (exists("GSEA") & !is.null(GSEA)) {
    if (dim(GSEA)[1] < 10){
    print (GSEA[,c("geneSet","description","normalizedEnrichmentScore","pValue","FDR")])
    } else {
      print(sprintf("There are %s significant Gene sets. Displaying only the %s most significant ones",dim(GSEA)[1],10))
      print (GSEA[1:10,c("geneSet","description","normalizedEnrichmentScore","pValue","FDR")])
    }
    rm(GSEA)
  } else{
    print ("No significant results...")
  }
}





