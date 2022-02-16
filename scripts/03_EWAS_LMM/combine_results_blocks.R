full_df=data.frame()
for (filename in snakemake@input[["results_blocks"]]){
  load(filename)
  full_df=rbind(full_df,results_df)
}

save(full_df,file=snakemake@output[["diff_meth_df"]])

