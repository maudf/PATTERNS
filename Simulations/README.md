# Run polygenic adaptation simulations using SLiM
Author: Maud Fagny

## Introduction
This document describe all the data files and scripts necessary to replicates figures, supplementary figures, supplementary tables and datasets from the following paper:  
Fagny M, Paulson JN, Kuijjer ML, Sonawane AR, Chen C.-Y., Lopes-Ramos CM, Glass K, Quackenbush J, Platig J. (2017) Exploring regulation in tissues with eQTL networks. _PNAS_ __114(37)__:E7841-E7850. [doi:10.1073/pnas.1707375114]()  

## General settings

###List of necessary softwares
Running the following scripts requires to have the following softwares and packages installed:

Software               | Where to find it  
---------------------- | ---------------------------------------------------  
plink2 v1.90 or higher | [https://www.cog-genomics.org/plink2](https://www.cog-genomics.org/plink2)  
htslib v1.12 or higher | [http://www.htslib.org](http://www.htslib.org)  
SelectionHapStats      | [https://github.com/ngarud/SelectionHapStats/](https://github.com/ngarud/SelectionHapStats)  

### List of required data files
Running these scripts do not require any additional data files.


## Scripts pipeline
To reproduce results and figures from our paper, run the following pipeline.
```{r rsetup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, tidy = FALSE)
```

### Run Simulations
The scripts in this section allow to run the simulations using SLiM.

#### Run Simulations under Neutrality
This script runs the simulations under neutrality using SLiM, with two populations of constant size 10,000, a fixed mutation rate of 2.38x10-8, a fixed mutation rate of 15x10-5. The simulations generate 30 independent 100kb-long sequences, with one QTL in the middle of the first 20 sequences, with a global frequency of 0.05. 100 individuals per populations are then sampled at generation 10,000 and their genotype is stored in a VCF file.
Two versions of the pipeline are proposed. This pipeline require to run a lot of occurences of the same script in order to get enough replicates for each simulation to be able to compute statistical power later on. Using a computer cluster is thus an easy way to compute these simulations. I thus propose a version of the pipeline that can be run on a cluster using SLURM as a workload manager. Here I show an example with 100 replicates: 
```{bash qc, eval=FALSE}
mkdir logs/
sbatch --array=1-100 code/wrap_simulations_Neutral_SLURM.sh -s 0.78 -f 0.05 -d ./ --overwrite
```
Another version of this script can be used directly on any computer with bash and SLiM installed:
```{bash qc, eval=FALSE}
bash code/wrap_simulations_Neutral.sh -r 100 -s 0.78 -f 0.05 -d ./ --overwrite
```

#### Run Simulations under Stabilizing Selection
In this simulation, the phenotype that is the sum of the effect size of each QTL in the first 20 chromosomes are under Stabilizing selection in the 2 populations. The selection start after 1000 generations. Here I show again an example with 100 replicates to be used on a SLURM-run cluster: 
```{bash qc, eval=FALSE}
sbatch --array=1-100 code/wrap_simulations_Stabilizing_SLURM.sh -o 0 -s 0.78 -f 0.05 -d ./ --overwrite
```
And on on any computer with bash and SLiM installed:
```{bash qc, eval=FALSE}
bash code/wrap_simulations_Stabilizing.sh -r 100 -o 0 -s 0.78 -f 0.05 -d ./ --overwrite
```

#### Run Simulations under Directional Selection
In this simulation, the phenotype that is the sum of the effect size of each QTL in the first 20 chromosomes is under Directional selection in the 2 populations with 2 different optimum: 0 in population 1 and another value in population 2. The selection start after 1000 generations. Here I show again examples with 100 replicates to be used on a SLURM-run cluster, with optimum for population 2 of 2, 5 and 10: 
```{bash qc, eval=FALSE}
sbatch --array=1-100 code/wrap_simulations_Directional_SLURM.sh -o 2 -s 0.78 -f 0.05 -d ./ --overwrite
sbatch --array=1-100 code/wrap_simulations_Directional_SLURM.sh -o 5 -s 0.78 -f 0.05 -d ./ --overwrite
sbatch --array=1-100 code/wrap_simulations_Directional_SLURM.sh -o 10 -s 0.78 -f 0.05 -d ./ --overwrite
```
And on on any computer with bash and SLiM installed:
```{bash qc, eval=FALSE}
bash code/wrap_simulations_Directional.sh -r 100 -o 2 -s 0.78 -f 0.05 -d ./ --overwrite
bash code/wrap_simulations_Directional.sh -r 100 -o 5 -s 0.78 -f 0.05 -d ./ --overwrite
bash code/wrap_simulations_Directional.sh -r 100 -o 10 -s 0.78 -f 0.05 -d ./ --overwrite
```

### Reformat output files

#### Reformat the vcf to run PCadapt
PCadapt takes BED plink files in input. The reformatting can be done using plink1.9.
```{bash qc, eval=FALSE}
bash code/convert_VCF_BED.sh -s 0.78 -f 0.05 -d ./Results/Neutral/
bash code/convert_VCF_BED.sh -o 0 -s 0.78 -f 0.05 -d ./Results/Stabilizing/
bash code/convert_VCF_BED.sh -o 2 -s 0.78 -f 0.05 -d ./Results/Directional/
bash code/convert_VCF_BED.sh -o 5 -s 0.78 -f 0.05 -d ./Results/Directional/
```

#### Reformat the vcf to run SelectionHapStats
SelectionHapStats takes a specific file format in input: comma-separated file with one snp on each line. The first column contains the SNP position. For phased VCF, the following columns contain the SNP allele for each haplotype (A, C, G, T, N). The reformatting was done using a bash routine. An example is given for each file type:
```{bash qc, eval=FALSE}
bash code/convert_VCF_BED_Neutral.sh -r 100 -s 0.78 -f 0.05 -d ./Results/Neutral/
bash code/convert_VCF_BED.sh -r 100 -o 0 -s 0.78 -f 0.05 -d ./Results/Stabilizing/
bash code/convert_VCF_BED.sh -r 100 -o 2 -s 0.78 -f 0.05 -d ./Results/Directional/
bash code/convert_VCF_BED.sh -r 100 -o 5 -s 0.78 -f 0.05 -d ./Results/Directional/
```

