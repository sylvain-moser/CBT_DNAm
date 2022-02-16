# Introduction
This is a snakemake pipeline developped to analyse changes in DNA-methylation during the course of Cognitive-Behavioural-Therapy with exposure intervention in Panic disorder patients.
This pipeline was used to perfrom all analysis described in Moser S, Martins J, Czamara D, Lange J, Müller-Myhsok B, Erhardt A. DNA-methylation dynamics across short-term, exposure-containing CBT in patients with panic disorder. Transl Psychiatry. 2022 Feb 1;12(1):46. doi: 10.1038/s41398-022-01802-7. PMID: 35105872.

# Installation: 
Snakemake runs the pipeline inside conda environements which it will install itself. The only requirements to run the pipeline is to have the conda environment management system installed : see https://docs.conda.io/projects/conda/en/latest/index.html

## Analysis on full dataset: 
The Epigenome-wide analysis is set-up to be run on a cluster with the SLURM job submission system using the following code: 

1) First, create the PD_exposure_main environment if not already created with:

```
conda env create -f envs/PD_exposure_main_env.yml
```

2) activate the PD_exposure_main environment and run the analysis 
```
conda activate PD_exposure_main
./run_analysis.sh
```

It produces all the figures and tables shown in Moser and al. It necessitates the actual patient data which cannot be shared online, but can be shared upon reasonable request to the authors.

## Example data-set analysis: 
For illustration and testing purposes, the analysis can be run with a smaller and randomly simulated dataset available in this repo. This anlalysis is set-up to be run on a personal computer and should take ca 1-2h. It produces the same type of outpout as the analysis on the real dataset, but the data being random, the results are biologically not meaningfull. 

The analysis can be run locally on a personal computer with the following code: 
1) First, create the PD_exposure_main environment if not already created with:

```
conda env create -f envs/PD_exposure_main_env.yml
```

2) activate the PD_exposure_main environment and run the analysis 
```
conda activate PD_exposure_main
./run_example_local.sh
```

It can be run on a cluster with the SLURM job submission system using the following code: 

1) First, create the PD_exposure_main environment if not already created with:

```
conda env create -f envs/PD_exposure_main_env.yml
```

2) activate the PD_exposure_main environment and run the analysis 
```
conda activate PD_exposure_main
./run_example.sh
```