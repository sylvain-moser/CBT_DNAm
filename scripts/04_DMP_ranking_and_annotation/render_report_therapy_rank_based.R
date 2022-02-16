library(rmarkdown)

rmarkdown::render(as.character(snakemake@params[["script_path"]]),params = list(
  therapy_final_processed_results=paste0("../../",snakemake@input[["therapy_final_processed_results"]]),
  therapy_cpg_ranked_anotated=paste0("../../",snakemake@input[["therapy_cpg_ranked_anotated"]]),
  WebGestal_script=paste0("../../",snakemake@params[["WebGestal_script"]]),
  back_to_normal_script=paste0("../../",snakemake@params[["back_to_normal_script"]]),
  methylation_df=paste0("../../",snakemake@input[["methylation_df"]]),
  methylation_residuals_df=paste0("../../",snakemake@input[["methylation_residuals_df"]]),
  therapy_cpg_ranked=paste0("../../",snakemake@input[["therapy_cpg_ranked"]]),
  candidates=paste0("../../",snakemake@input[["candidates"]])
),output_file =paste0("../../",snakemake@output[["report"]]))
