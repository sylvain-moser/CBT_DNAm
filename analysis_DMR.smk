############################################################################
# main script to run the DMR analysis on biological correlates of CBT          #
############################################################################
from os import path
import numpy as np

prefix_statgen = "/Users/sylvain_moser/psycl_statgen/symo/PD_redo" if path.exists("/Users/sylvain_moser/psycl_statgen/symo/dummy_file.txt")  else "/psycl/g/mpsstatgen/symo/PD_redo"

rule all:
    input:
        DMR_result_file=expand(os.path.join(prefix_statgen,"results","05_DMR","DMR_{condition}.regions-t.bed"),condition=["exposure","therapy"]),
        annotated_DMR=expand(os.path.join(prefix_statgen,"results","05_DMR","annotated_DMR_{condition}.txt"),condition=["exposure"])


rule prepare_combp_input:
    input:
        diff_meth_df=os.path.join(prefix_statgen,"results","04_DMP_ranking_and_annotation","regulated_island_{condition}_INT_processed.RData"),
        illumina_annot_file=os.path.join(prefix_statgen,"data","prepared_data","illumina_450k_annotations.RDS")
    output:
        regulated_island_bim_file=os.path.join(prefix_statgen,"results","05_DMR","DMR_{condition,[A-Za-z0-9]+}.bed")
    params:
        pval_column="best_lmm_pval" # a bit anti-conservative maybe? 
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join(prefix_statgen,"workflow","scripts","05_DMR","prepare_bed_with_pval.R")

rule find_DMR:
    input:
        regulated_island_bim_file=os.path.join(prefix_statgen,"results","05_DMR","DMR_{condition}.bed")
    output:
        DMR_result_file=os.path.join(prefix_statgen,"results","05_DMR","DMR_{condition,[A-Za-z0-9]+}.regions-t.bed")
    params:
        DMR_file_prefix=os.path.join(prefix_statgen,"results","05_DMR","DMR_{condition}")
    conda:
        "envs/PD_exposure_combp_env.yml"
    shell:
        "comb-p pipeline -c 4 --seed 0.05 --dist 750 -p {params.DMR_file_prefix} --region-filter-p 0.1 {input.regulated_island_bim_file}" #best parameters reported in Mallik 2019

rule annotate_DMR:
    input:
        DMR_result_file=os.path.join(prefix_statgen,"results","05_DMR","DMR_{condition}.regions-t.bed")
    params:
        n_probe_threshold=2
    output:
        DMR_near_promoter=os.path.join(prefix_statgen,"results","05_DMR","genes_near_promoter_{condition,[A-Za-z0-9]+}.txt"),
        annotated_DMR_near_promoter=os.path.join(prefix_statgen,"results","05_DMR","annotated_DMR_near_promoter_{condition}.txt"),
        DMR=os.path.join(prefix_statgen,"results","05_DMR","genes_{condition}.txt"),
        annotated_DMR=os.path.join(prefix_statgen,"results","05_DMR","annotated_DMR_{condition}.txt")
    conda:
        "envs/PD_exposure_r4_env.yml"
    script:
        os.path.join(prefix_statgen,"workflow","scripts","05_DMR","annotate_DMR.R")

