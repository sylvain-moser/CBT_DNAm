library(rmarkdown)

rmarkdown::render(as.character(snakemake@params[["script_path"]]),params = list(
  exposure_final_processed_results=paste0("../../",snakemake@input[["exposure_final_processed_results"]]),
  exposure_cpg_ranked_anotated=paste0("../../",snakemake@input[["exposure_cpg_ranked_anotated"]]),
  WebGestal_script=paste0("../../",snakemake@params[["WebGestal_script"]]),
  back_to_normal_script=paste0("../../",snakemake@params[["back_to_normal_script"]]),
  methylation_df=paste0("../../",snakemake@input[["methylation_df"]]),
  methylation_residuals_df=paste0("../../",snakemake@input[["methylation_residuals_df"]]),
  exposure_cpg_ranked=paste0("../../",snakemake@input[["exposure_cpg_ranked"]]),
  candidates=paste0("../../",snakemake@input[["candidates"]])
),output_file =paste0("../../",snakemake@output[["report"]]))
