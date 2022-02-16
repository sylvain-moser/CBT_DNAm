############################################################################
# main script to run the analysis on biological correlates of CBT          #
############################################################################
from os import path
import numpy as np

##### original basis data file ######

NID_timepoints_keys_file_path=os.path.join("data_example","initial_data","M00832_Platte1-3_samplesheet_example_random.sav") 
cell_count_file_path=os.path.join("data_example","initial_data","CellTypeComposition_random.txt")
pheno_file_path=os.path.join("data_example","initial_data","pheno_example_random.csv")
methylation_file_path=os.path.join("data_example","initial_data","filtered_meth_random.RData")
expression_file_path=os.path.join("data_example","initial_data","gene_expression_example_random.RData")
Sample_ID_keys_path=os.path.join("data_example","initial_data","sampleSheet_pheno_example_random.txt")
cell_composition_expression_file_path=os.path.join("data_example","initial_data","CellComp_SVPS_random.tab")


rule all:
    input:
        cell_types_results_df=expand(os.path.join("results_example","02_immune_cells","{condition}_LMMs_results.RDS"),condition=["exposure","therapy"]),
        CD4T_plot=expand(os.path.join("results_example","02_immune_cells/CD4T_cells_{condition}.png"),condition=["exposure","therapy"]),
        QQplots=expand(os.path.join("results_example","03_EWAS_LMM","QQplots_{condition}.png"),condition=["exposure","therapy"]),
        general_report=expand(os.path.join("results_example","04_DMP_ranking_and_annotation","{condition}_rank_based_report.html"),condition=["exposure","therapy"]),
        candidate_analysis=os.path.join("results_example","06_candiate_gene_analysis","candidate_analysis.html"),
        meth_profiles=expand(os.path.join("results_example","07_methylation_profiles","{gene}_CBT_methylation_profile.html"),gene=["HTR3A"]),
        expression_profiles=expand(os.path.join("results_example","08_expression_profiles","{gene}_CBT_expression_profile.html"),gene=["HTR3A"]), #"KLF16","TMEM132D"
        cell_types_results_plots=os.path.join("results_example","02_immune_cells","therapy_LMM_T0_E_plots.png"),
        figure2=os.path.join("results_example","02_immune_cells","figure2.pdf"),
        figure3=os.path.join("results_example","04_DMP_ranking_and_annotation","figure3.eps")


###### prepare the different data files for later analysis ######

rule prepare_cell_types_df:
    input:
        NID_timepoints_keys_file=NID_timepoints_keys_file_path,
        cell_count_file=cell_count_file_path,
        pheno_file=pheno_file_path
    output:
        cell_types_df=os.path.join("data_example","prepared_data","cell_types_df.RDS")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts","01_prepare_data","prepare_cell_types_df.R")

rule prepare_meth_df:
    input:
        methylation_file=methylation_file_path,
        NID_timepoints_keys_file=NID_timepoints_keys_file_path,
        cell_count_file=cell_count_file_path,
        pheno_file=pheno_file_path,
    output:
        methylation_df=os.path.join("data_example","prepared_data","methylation_df_cellcount_smoking.Rdata"),
        methylation_df_feather=os.path.join("data_example","prepared_data","methylation_df_cellcount_smoking.feather")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts/01_prepare_data/Match_Sentrix_ID_to_NID_timepoints.R")

rule compute_methylation_residuals: 
    input:
        methylation_df=os.path.join("data_example","prepared_data","methylation_df_cellcount_smoking.Rdata")
    output:
        methylation_residuals_df=os.path.join("data_example","prepared_data","methylation_residuals.Rdata")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts/01_prepare_data/compute_methylation_residuals.R")

rule prepare_expression_df:
    input:
        expression_file=expression_file_path,
        Sample_ID_keys=Sample_ID_keys_path,
        NID_timepoints_keys_file=NID_timepoints_keys_file_path,
        cell_composition_expression=cell_composition_expression_file_path,
        pheno_file=pheno_file_path
    output:
        expression_df=os.path.join("data_example","prepared_data","gene_expression_df.RDS")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts/01_prepare_data/prepare_expression_df.R")

rule compute_expression_residuals: 
    input: 
        expression_df=os.path.join("data_example","prepared_data","gene_expression_df.RDS")
    output:
        expression_residuals_df=os.path.join("data_example","prepared_data","expression_residuals.RDS")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts/01_prepare_data/compute_expression_residuals.R")

rule prepare_illumina_annotations:
    input: 
        []
    conda: 
        "envs/PD_exposure_r_illumina450k_env.yml"
    output:
        illumina_annot_file=os.path.join("data_example","prepared_data","illumina_450k_annotations.RDS")
    script:
        os.path.join("scripts/01_prepare_data/make_illumina_annotations_df.R")

################## Immune cell-types Analysis ##################

rule compute_cell_types_LMM_exposure:
    input:
        cell_types_df=os.path.join("data_example","prepared_data","cell_types_df.RDS"),
        pheno_file=pheno_file_path
    output:
        cell_types_results_df=os.path.join("results_example","02_immune_cells","exposure_LMMs_results.RDS"),
        cell_types_results_table=os.path.join("results_example","02_immune_cells","exposure_LMMs_results.txt"),
        immune_cells_features_df=os.path.join("results_example","02_immune_cells","immune_cells_exposure_features.csv")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts","02_immune_cells","LMM_immune_cell_counts_exposure.R")


rule compute_cell_types_LMM_therapy:
    input:
        cell_types_df=os.path.join("data_example","prepared_data","cell_types_df.RDS"),
        pheno_file=pheno_file_path
    output:
        cell_types_results_df=os.path.join("results_example","02_immune_cells","therapy_LMMs_results.RDS"),
        cell_types_results_table=os.path.join("results_example","02_immune_cells","therapy_LMMs_results.txt"),
        immune_cells_features_df=os.path.join("results_example","02_immune_cells","immune_cells_therapy_features.csv")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts","02_immune_cells","LMM_immune_cell_counts_therapy.R")


rule compute_cell_types_LMM_T0_E_therapy:
    input:
        cell_types_df=os.path.join("data_example","prepared_data","cell_types_df.RDS"),
        pheno_file=pheno_file_path
    output:
        cell_types_results_table=os.path.join("results_example","02_immune_cells","therapy_LMM_T0_E_results.txt"),
        cell_types_results_plots=os.path.join("results_example","02_immune_cells","therapy_LMM_T0_E_plots.png"),
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts","02_immune_cells","LMM_T0_E_immune_cells_proportions_therapy.R")

rule plot_cell_types:
    input:
        cell_types_df=os.path.join("data_example","prepared_data","cell_types_df.RDS"),
        pheno_file=pheno_file_path
    output:
        CD4T_plot=os.path.join("results_example","02_immune_cells/CD4T_cells_{condition}.png")
    params:
        plot_dir=os.path.join("results_example","02_immune_cells")
    conda:
        "envs/PD_exposure_r_env.yml"
    script: 
        os.path.join("scripts","02_immune_cells","plot_cells_proportions_{wildcards.condition}.R")

rule make_figure2:
    input:
        cell_types_df=os.path.join("data_example","prepared_data","cell_types_df.RDS"),
        pheno_file=pheno_file_path
    output:
        figure2_path=os.path.join("results_example","02_immune_cells","figure2.pdf")
    conda:
        "envs/PD_exposure_r_figures_env.yml"
    script:
        os.path.join("scripts","02_immune_cells","figure2_boxplot_and_line.R")

################# Methylome-wide search of regulated CpGs##########

# Step1: Computation of the LMMs 
rule find_DMP_INT:
    input:
        methylation_df=os.path.join("data_example","prepared_data","methylation_df_cellcount_smoking.Rdata")
    params:
        start_idx="{start_idx}",
        surrogate_variable_script=os.path.join("scripts","03_EWAS_LMM","compute_surrogate_variables.R"),
        surrogate_variable_matrix=os.path.join("results_example","03_EWAS_LMM","surrogate_matrix_{condition}.RData"),
        surrogate_df_path=os.path.join("results_example","03_EWAS_LMM","surrogate_variables_{condition}/"),
        output_dir=os.path.join("results_example","03_EWAS_LMM","regulated_island_{condition}_")
    output:
        find_DMP_results=os.path.join("results_example","03_EWAS_LMM","regulated_island_{condition}_{start_idx}_second_lmertestREMLINT.Rdata"),
        dummy_output=os.path.join("results_example","03_EWAS_LMM","surrogate_variables_{condition}","dummy_{condition}_{start_idx}")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts","03_EWAS_LMM","find_differentially_methylated_random_slope_{wildcards.condition}_lmertestREMLINT_example.R")

# Step2: Building Diagnositcs plots

rule make_LMM_qqplots: 
    input:
        LMM_results=expand(os.path.join("results_example","03_EWAS_LMM","regulated_island_{{condition}}_{chunk}_second_lmertestREMLINT.Rdata"),chunk=range(1))
    output:
        QQplots=os.path.join("results_example","03_EWAS_LMM","QQplots_{condition}.png")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts/03_EWAS_LMM/make_qqplots_exposure_lmerREMLINT_figure.R")

rule combine_diff_regulated_block_results:
    input:
        results_blocks=expand(os.path.join("results_example","03_EWAS_LMM","regulated_island_{{condition}}_{block}_second_lmertestREMLINT.Rdata"),block=range(1))
    output:
        diff_meth_df=os.path.join("results_example","03_EWAS_LMM","regulated_island_{condition}_INT_full.RData")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts","03_EWAS_LMM","combine_results_blocks.R")

# #Step 3:  Ranking and And annotation of the Best CpGs

rule process_island_results:
    input: 
        diff_meth_df=os.path.join("results_example","03_EWAS_LMM","regulated_island_{condition}_INT_full.RData")
    output:
        processed_diff_meth_df=os.path.join("results_example","04_DMP_ranking_and_annotation","regulated_island_{condition}_INT_processed.RData"),
        processed_diff_meth_df_feather=os.path.join("results_example","04_DMP_ranking_and_annotation","regulated_island_{condition}_INT_processed.feather"),
        failed_models=os.path.join("results_example","04_DMP_ranking_and_annotation","failed_islands_{condition}_INT.RData")
    conda:
        "envs/PD_exposure_r_env.yml"
    script: 
        os.path.join("scripts","04_DMP_ranking_and_annotation","process_results_dplyr.R")

rule rank_cpgs:
    input:
        processed_diff_meth_df=os.path.join("results_example","04_DMP_ranking_and_annotation","regulated_island_{condition}_INT_processed.RData"),
        meth_df=os.path.join("data_example","prepared_data","methylation_df_cellcount_smoking.Rdata"),
        illumina_annot_file=os.path.join("data_example","prepared_data","illumina_450k_annotations.RDS")
    output:
        ranked_cpgs=os.path.join("results_example","04_DMP_ranking_and_annotation","{condition}_INT_ranked.RData"),
        ranked_annotated_cpgs=os.path.join("results_example","04_DMP_ranking_and_annotation","{condition}_INT_ranked_annotated.RData"),
        ranked_annotated_cpgs_csv=os.path.join("results_example","04_DMP_ranking_and_annotation","{condition}_INT_ranked_annotated.csv")
    params:
        cores=20
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts","04_DMP_ranking_and_annotation","compute_ranking_parallel.R")

rule render_report_exposure:
    input:
        exposure_final_processed_results=os.path.join("results_example","04_DMP_ranking_and_annotation","regulated_island_exposure_INT_processed.RData"),
        exposure_cpg_ranked_anotated=os.path.join("results_example","04_DMP_ranking_and_annotation","exposure_INT_ranked_annotated.RData"),
        methylation_df=os.path.join("data_example","prepared_data","methylation_df_cellcount_smoking.Rdata"),
        methylation_residuals_df=os.path.join("data_example","prepared_data","methylation_residuals.Rdata"),
        exposure_cpg_ranked=os.path.join("results_example","04_DMP_ranking_and_annotation","exposure_INT_ranked.RData"),
        candidates=os.path.join("scripts","06_candiate_gene_analysis","candidates.txt")
    params:
        WebGestal_script=os.path.join("scripts","04_DMP_ranking_and_annotation","webgestalt_GSEA_450k.R"),
        back_to_normal_script=os.path.join("scripts","04_DMP_ranking_and_annotation","find_cpgs_not_returning_to_norma_exposure.R"),
        script_path=os.path.join("scripts","04_DMP_ranking_and_annotation","final_report_exposure_rank_based.Rmd")
    output:
        report=os.path.join("results_example","04_DMP_ranking_and_annotation","exposure_rank_based_report.html")
    conda:
        "envs/PD_exposure_r4_env.yml"
    script:
        os.path.join("scripts","04_DMP_ranking_and_annotation","render_report_exposure_rank_based.R")

rule render_report_therapy:
    input:
        therapy_final_processed_results=os.path.join("results_example","04_DMP_ranking_and_annotation","regulated_island_therapy_INT_processed.RData"),
        therapy_cpg_ranked_anotated=os.path.join("results_example","04_DMP_ranking_and_annotation","therapy_INT_ranked_annotated.RData"),
        methylation_df=os.path.join("data_example","prepared_data","methylation_df_cellcount_smoking.Rdata"),
        methylation_residuals_df=os.path.join("data_example","prepared_data","methylation_residuals.Rdata"),
        therapy_cpg_ranked=os.path.join("results_example","04_DMP_ranking_and_annotation","therapy_INT_ranked.RData"),
        candidates=os.path.join("scripts","06_candiate_gene_analysis","candidates.txt")
    params:
        WebGestal_script=os.path.join("scripts","04_DMP_ranking_and_annotation","webgestalt_GSEA_450k.R"),
        back_to_normal_script=os.path.join("scripts","04_DMP_ranking_and_annotation","find_cpgs_not_returning_to_norma_therapy.R"),
        script_path=os.path.join("scripts","04_DMP_ranking_and_annotation","final_report_therapy_rank_based.Rmd")
    output:
        report=os.path.join("results_example","04_DMP_ranking_and_annotation","therapy_rank_based_report.html")
    conda:
        "envs/PD_exposure_r4_env.yml"
    script:
        os.path.join("scripts","04_DMP_ranking_and_annotation","render_report_therapy_rank_based.R")

rule make_figure3:
    input:
        exposure_cpg_ranked_anotated=os.path.join("results_example","04_DMP_ranking_and_annotation","exposure_INT_ranked_annotated.RData"),
        exposure_cpg_ranked=os.path.join("results_example","04_DMP_ranking_and_annotation","exposure_INT_ranked.RData"),
        therapy_cpg_ranked_anotated=os.path.join("results_example","04_DMP_ranking_and_annotation","therapy_INT_ranked_annotated.RData"),
        therapy_cpg_ranked=os.path.join("results_example","04_DMP_ranking_and_annotation","therapy_INT_ranked.RData"),
        methylation_residuals_df=os.path.join("data_example","prepared_data","methylation_residuals.Rdata"),
        pheno_file=pheno_file_path,
        expression_residuals_df=os.path.join("data_example","prepared_data","expression_residuals.RDS")
    output:
        figure3_path=os.path.join("results_example","04_DMP_ranking_and_annotation","figure3.eps")
    conda:
        "envs/PD_exposure_r_figures_env.yml"
    script:
        os.path.join("scripts","04_DMP_ranking_and_annotation","figure3.R")

# Gene candidate Analysis: 

rule candidate_analysis:
    input: 
        candidates=os.path.join("scripts","06_candiate_gene_analysis","candidates.txt"),
        diff_meth_df=os.path.join("results_example","03_EWAS_LMM","regulated_island_{condition}_INT_full.RData"),
        illumina_annot_file=os.path.join("data_example","prepared_data","illumina_450k_annotations.RDS")
    output:
        candidate_analysis_results=os.path.join("results_example","06_candiate_gene_analysis","candidate_genes_analysis_{condition,[A-Za-z]+}_results.txt")
    conda:
        "envs/PD_exposure_r4_env.yml"
    script:
        os.path.join("scripts","06_candiate_gene_analysis","gene_candidate_analysis.R")


rule candidate_analysis_expression:
    input: 
        candidates=os.path.join("scripts","06_candiate_gene_analysis","candidates_expression.csv"),
        diff_meth_df=os.path.join("results_example","03_EWAS_LMM","regulated_island_{condition}_INT_full.RData"),
        illumina_annot_file=os.path.join("data_example","prepared_data","illumina_450k_annotations.RDS")
    output:
        candidate_analysis_results=os.path.join("results_example","06_candiate_gene_analysis","candidate_genes_analysis_expression_{condition}_results.txt")
    conda:
        "envs/PD_exposure_r4_env.yml"
    script:
        os.path.join("scripts","06_candiate_gene_analysis","gene_candidate_analysis_expression.R")

rule candidate_analysis_iurato_cpgs:
    input: 
        candidates=os.path.join("scripts","06_candiate_gene_analysis","candidates.txt"),
        diff_meth_df=os.path.join("results_example","03_EWAS_LMM","regulated_island_{condition}_INT_full.RData"),
        illumina_annot_file=os.path.join("data_example","prepared_data","illumina_450k_annotations.RDS")
    output:
        candidate_analysis_results=os.path.join("results_example","06_candiate_gene_analysis","candidate_genes_analysis_iurato_{condition,[A-Za-z]+}_results.txt")
    conda:
        "envs/PD_exposure_r4_env.yml"
    script:
        os.path.join("scripts","06_candiate_gene_analysis","gene_candidate_analysis_iurato_cpgs.R")

rule render_report_candidate_analysis:
    input: 
        rmd_script=os.path.join("scripts","06_candiate_gene_analysis","candidate_analysis.Rmd"),
        candidate_analysis_exposure_results=os.path.join("results_example","06_candiate_gene_analysis","candidate_genes_analysis_exposure_results.txt"),
        exposure_final_processed_results=os.path.join("results_example","04_DMP_ranking_and_annotation","regulated_island_exposure_INT_processed.RData"),
        candidate_analysis_therapy_results=os.path.join("results_example","06_candiate_gene_analysis","candidate_genes_analysis_therapy_results.txt"),
        therapy_final_processed_results=os.path.join("results_example","04_DMP_ranking_and_annotation","regulated_island_therapy_INT_processed.RData"),
        candidate_analysis_expression_therapy_results=os.path.join("results_example","06_candiate_gene_analysis","candidate_genes_analysis_expression_therapy_results.txt"),
        candidate_analysis_expression_exposure_results=os.path.join("results_example","06_candiate_gene_analysis","candidate_genes_analysis_expression_exposure_results.txt"),
        candidate_analysis_iurato_exposure_results=os.path.join("results_example","06_candiate_gene_analysis","candidate_genes_analysis_iurato_exposure_results.txt"),
        candidate_analysis_iurato_therapy_results=os.path.join("results_example","06_candiate_gene_analysis","candidate_genes_analysis_iurato_therapy_results.txt")
    output:
        report=os.path.join("results_example","06_candiate_gene_analysis","candidate_analysis.html")
    conda:
        "envs/PD_exposure_r4_env.yml"
    script:
        os.path.join("scripts","06_candiate_gene_analysis","render_report_candidate_analysis.R")
    

# Step 4: Methylation Profile of the Cpgs of interest: 

rule render_meth_profile:
    input:
        methylation_df=os.path.join("data_example","prepared_data","methylation_df_cellcount_smoking.Rdata"),
        pheno=pheno_file_path,
        methylation_residuals=os.path.join("data_example","prepared_data","methylation_residuals.Rdata"),
        methylation_profile_script=os.path.join("scripts","07_methylation_profiles","{gene}_CBT_profile.Rmd")
    output:
        report=os.path.join("results_example","07_methylation_profiles","{gene}_CBT_methylation_profile.html")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts","07_methylation_profiles","render_methylation_profile.R")


# Step 5: Gene expression Profile of the CpGs of interest:

rule render_expression_profile:
    input:
        expression_residuals_df=os.path.join("data_example","prepared_data","expression_residuals.RDS"),
        pheno=pheno_file_path,
        illu_annot=os.path.join("data_example","initial_data","HumanHT-12_V4_0_R2_15002873_B_example.txt"),
        expression_profile_script=os.path.join("scripts","08_expression_profiles","{gene}_expression_CBT_profil.Rmd")
    output:
        report=os.path.join("results_example","08_expression_profiles","{gene}_CBT_expression_profile.html")
    conda:
        "envs/PD_exposure_r_env.yml"
    script:
        os.path.join("scripts","08_expression_profiles","render_expression_profile.R")

