library(rmarkdown)

rmarkdown::render(as.character(snakemake@input[["rmd_script"]]),params = list(
  candidate_analysis_exposure_results=paste0("../../",snakemake@input[["candidate_analysis_exposure_results"]]),
  exposure_final_processed_results=paste0("../../",snakemake@input[["exposure_final_processed_results"]]),
  candidate_analysis_therapy_results=paste0("../../",snakemake@input[["candidate_analysis_therapy_results"]]),
  therapy_final_processed_results=paste0("../../",snakemake@input[["therapy_final_processed_results"]]),
  candidate_analysis_expression_therapy_results=paste0("../../",snakemake@input[["candidate_analysis_expression_therapy_results"]]),
  candidate_analysis_expression_exposure_results=paste0("../../",snakemake@input[["candidate_analysis_expression_exposure_results"]]),
  candidate_analysis_iurato_therapy_results=paste0("../../",snakemake@input[["candidate_analysis_iurato_therapy_results"]]),
  candidate_analysis_iurato_exposure_results=paste0("../../",snakemake@input[["candidate_analysis_iurato_exposure_results"]])
),output_file =paste0("../../",snakemake@output[["report"]]))
