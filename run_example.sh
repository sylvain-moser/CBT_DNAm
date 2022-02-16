mkdir -p snakemake_logs_example ; snakemake --snakefile analysis_example.smk all -j 10 --use-conda --cluster-config cluster_example.json --cluster "sbatch --time={cluster.time} --mem={cluster.mem} --exclude={cluster.excl} -o {cluster.stdout} -e {cluster.stderr} --job-name={cluster.jobname} --cpus-per-task={cluster.cpus}" > run_example.out 2> run_example.err