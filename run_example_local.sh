mkdir -p snakemake_logs_example ;snakemake --snakefile analysis_example.smk all --use-conda --cores 4 > run_example.out 2> run_example.err
