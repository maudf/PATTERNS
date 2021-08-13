#!/usr/bin/Rscript
### Maud Fagny
### 2021-08-13
### makeCondorMatrixPlot.R
### Make a matrix representation of each eQTL network
###-----------------------------------------------------


### Load libraries
require(condor)
require(igraph)
require(getopt)
require(data.table)
require(this.path)

### Set variables
spec = matrix(c('help','h', 0, "logical",
                'directory', 'd', 1, "character",
                'output', 'o', 1, "character"),
              byrow=TRUE, ncol=4)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

if ( is.null(opt$directory    ) ) {opt$directory    = './Results/Networks/'     }
if ( is.null(opt$output    ) ) {opt$output    = './Results/GO.results/'     }

cat("Running compute_stats.R with options:\n",
    "directory =", opt$directory, "\n",
    "output =", opt$output, "\n")

### Load necessary functions 
script.path=gsub("/[^/]*$", "/", this.path())
source(paste0(script.path, 'condor.GO.R'))

### Create output folder if not exists
if (!(dir.exists(opt$output))){dir.create(opt$output, recursive = T, showWarnings = FALSE)}

### Enrichment in GO terms among communities
go.result <- 
tab.GO <- NULL
files=list.files(path=opt$directory, pattern = '.*.rds')
for (f in paste0(opt$directory, files)){
    tissue=gsub(".*/Edges_matrixeqtl_", "", gsub("_clusters.rds", "", f))
    cat("Load network for", tissue, "\n")
    condor.modularity=readRDS(f)

    cat("Computing GO and KEGG enrichment analyses for", tissue, "...\n")
    condor.modularity$qscores$blue.qscore$blue.names <- gsub("\\.[0-9][0-9]*", "",
                                                             condor.modularity$qscores$blue.qscore$blue.names)
    condor.modularity$qscores$red.qscore$red.names <- gsub("\\.[0-9][0-9]*", "",
                                                           condor.modularity$qscores$red.qscore$red.names)
    condor.modularity$blue.memb$blue.names <-gsub("\\.[0-9][0-9]*", "",
                                                           condor.modularity$blue.memb$blue.names)
    
    all.gene <- condor.modularity$qscores$blue.qscore$blue.names

    go.result <- condor.GO(condor.object = condor.modularity,
                                       go.ontology = c('BP', 'MF', 'CC'), unadj.p.cut = 0.1,
                                       min.overlap=5, symbol.universe = all.gene)
    rm(condor.modularity)
    saveRDS(go.result, file=paste0(opt$output, "GO.results.", tissue, ".rds"))
    tab.GO <- rbind(tab.GO, 
                    data.frame("Tissue"=rep(ti, nrow(go.result)), go.result, stringsAsFactor=F))
    
}
### Write results in a TAB-separated file
write.table(tab.GO, paste0( opt$output, "GO.results.", tissue, ".txt"),
            row.names=F, sep="\t", quote=F)
