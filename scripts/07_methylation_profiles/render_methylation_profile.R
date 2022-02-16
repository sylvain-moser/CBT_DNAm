library(rmarkdown)
rmarkdown::render(as.character(snakemake@input["methylation_profile_script"]),params = list(
  methylation_df = paste0("../../",snakemake@input[["methylation_df"]]),
  pheno =paste0("../../",snakemake@input[["pheno"]]),
  methylation_residuals=paste0("../../",snakemake@input[["methylation_residuals"]])
),output_file =paste0("../../",snakemake@output[["report"]]))
