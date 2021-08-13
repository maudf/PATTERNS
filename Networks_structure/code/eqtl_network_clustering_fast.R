#!/usr/bin/Rscript
### Maud Fagny
### 2021-08-13
### eqtl_network_clustering_fast.R
### Find module structure within eQTL networks
###-----------------------------------------------------


### Load libraries
require(condor)
require(getopt)
require(igraph)

### Set variables
spec = matrix(c('help','h', 0, "logical",
                'directory','d', 1, "character",
                'threshold','t', 1, "numeric",
                'output','o', 1, "character"),
              byrow=TRUE, ncol=4)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)}

if ( is.null(opt$directory    ) ) {opt$directory    = './eqtls/'     }
if ( is.null(opt$output    ) ) {opt$output    = './Networks/'     }
if ( is.null(opt$threshold    ) ) {opt$threshold    = 0.05     }

cat("Running eqtl_network_clustering_fast.R with options:\n",
    "directory =", opt$directory, "\n",
    "output =", opt$output, "\n",
    "threshold=", opt$threshold, "\n")

### Check if input folder exists
if(!dir.exists(opt$directory)) {stop(paste(opt$directory, "not found. Exit...\n"))}

### Create output folder if do not exists
if(!dir.exists(opt$output)){dir.create(opt$output, recursive = T)}

### List eqtls files

eqtl.files <- list.files(path=opt$directory, pattern = "*.[Rr][Dd][Ss]") # File containing eQTLs
if(length(eqtl.files)==0) {stop(opt$directory, "do not contains RDS files, Exit...\n")}

### Create eqtl network clusters 
for (f in eqtl.files){
  cat("___________________________________________________\n")
  cat("Reading", f, "...\n")
  eqtls <- readRDS(paste0(opt$directory,f))
  eqtls <- eqtls[eqtls$FDR<=opt$threshold,]
  
  cat("Preparing data for clustering analysis...\n")
  basename <- gsub(".[Rr][Dd][Ss]", "", f)
  elist<- data.frame("red" = eqtls$snps, "blue" = eqtls$gene)
  condor.object <- create.condor.object(elist)
  
  cat("Clustering network...\n")
  network <- condor.cluster(condor.object)
  output.file=paste0(opt$output, basename, "_clusters.rds")
  cat("Saving network in", output.file, "...\n")
  save(network, file=output.file)
  
  cat("Analysing modularity...\n")
  modularity <- condor.qscore(network)
  modularity.file=paste0(opt$output, basename, "_clusters.rds")
  cat("Saving modularity results in", modularity.file, "...\n")
  saveRDS(modularity, file=modularity.file)
}




