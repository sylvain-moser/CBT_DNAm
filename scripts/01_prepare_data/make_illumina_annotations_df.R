library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

data("Locations")
bed_core=as_tibble(Locations,rownames="island") %>% mutate(stop=pos+1) %>% dplyr::rename(start=pos) 
saveRDS(bed_core,file=snakemake@output[["illumina_annot_file"]])