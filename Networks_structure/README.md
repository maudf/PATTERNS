# Build and analyse structure of eQTL networks

## Introduction
This document describes all the data files and scripts necessary to
build and analyse the structures of eQTL networks as done in the PATTERNS project.

## General settings

### List of necessary R packages and softwares
Running the following scripts requires to have the following softwares and packages installed:

#### Softwares
Software               | Where to find it 
---------------------- | --------------------------------------------------- 
R                      | [https://cran.r-project.org/](https://cran.r-project.org/) 

#### R packages
Package        | Where to find it 
-------------- | ---------------------------------------------------------- 
getopt         | [https://cran.r-project.org/web/packages/getopt/index.html](https://cran.r-project.org/web/packages/getopt/index.html) 
RColorBrewer   | [https://cran.r-project.org/web/packages/RColorBrewer/index.html](https://cran.r-project.org/web/packages/RColorBrewer/index.html) 
R Bioconductor | [http://bioconductor.org](http://bioconductor.org) 
condor         | [https://github.com/jplatig/condor] (https://github.com/jplatig/condor)
igraph          |  [https://cran.r-project.org/web/packages/igraph/index.html] (https://cran.r-project.org/web/packages/igraph/index.html)
data.table    |  [https://cran.r-project.org/web/packages/data.table/index.html] (https://cran.r-project.org/web/packages/data.table/index.html)
plyr          |  [https://cran.r-project.org/web/packages/plyr/index.html](https://cran.r-project.org/web/packages/plyr/index.html)
this.path          |  [https://cran.r-project.org/web/packages/this.path/index.html](https://cran.r-project.org/web/packages/this.path/index.html)

#### R Bioconductor packages
Package        | Where to find it 
-------------- | -----------------------------------------------------------  
Biobase        | [http://bioconductor.org/packages/release/bioc/html/Biobase.html](http://bioconductor.org/packages/release/bioc/html/Biobase.html)
GOstats       | [https://bioconductor.org/packages/release/bioc/html/GOstats.html](https://bioconductor.org/packages/release/bioc/html/GOstats.html)
org.Hs.eg.db | [https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

## required files
Running these scripts require to have the following data files as input:
eQTL results stored as an RDS file containing an R data frame with at
least 3 colums:
* snps: SNPs ID
* gene: Genes ID
* FDR: FDR-corrected p-values
All the eqtl files to be analysed must be stored in the same folder.

## Pipeline 
### Build eQTL Networks using the R condor Package
This bash script runs the genotype QC and generates the input files for matrix eQTL. 
You can fix the FDR threshold above which eQTLs will not be included in the analysis.
**Warning:** This scripts can be very expensive in terms of memory usage. For 3,235,571 potential cis- and trans-eQTLs, plan for at least 200Go memory.
```{bash networks}
Rscript code/eqtl_network_clustering_fast.R -d ./eqtls -o ./Networks/ 
```

## Modularity Analysis
This bash script runs an R script that extract the modularity for each network:
```{bash modularity}
code/modularity.R -d ./Networks/ 
```

## Plot Network as a matrix
This bash script runs an R script that plot each network as a matrix with SNPs on x-axis and Gene on y-axis:
```{bash modularity}
code/makeCondorMatrixPlot.R -d ./Networks/ -o ./Figures/
```

## Modules Gene Ontology Enrichment Analysis
This bash script run an R script that will compute Gene Ontology enrichment analysis for each module in each network:
```{bash GO}
Rscript code/analyze_clusters_GO.R -d ./Networks/ -o ./GO.results/ 
```
