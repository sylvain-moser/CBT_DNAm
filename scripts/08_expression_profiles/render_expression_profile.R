library(rmarkdown)
rmarkdown::render(as.character(snakemake@input["expression_profile_script"]),params = list(
  expression_df=paste0("../../",snakemake@input[["expression_residuals_df"]]),
  pheno =paste0("../../",snakemake@input[["pheno"]]),
  illu_annot=paste0("../../",snakemake@input[["illu_annot"]])
),output_file =paste0("../../",snakemake@output[["report"]]))
