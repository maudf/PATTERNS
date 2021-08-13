#!/usr/bin/Rscript
### Maud Fagny
### 2021-08-13
### modularity.R
### Extract modularity from network files
###-----------------------------------------------------

### Load packages
require("getopt")

### Set variables
spec = matrix(c('help','h', 0, "logical",
                'directory','d', 1, "character"),
              byrow=TRUE, ncol=4)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$directory    ) ) {opt$directory    = './Results/Networks/'     }

cat("Running compute_stats.R with options:\n",
    "directory =", opt$directory, "\n")

### Load networks and extract modularity
mod=c()
clusters=c()
files=list.files(path=opt$directory, pattern = '.*.rds')
for (f in paste0(opt$directory, files)){
  tissue=gsub(".*/Edges_matrixeqtl_", "", gsub("_clusters.rds", "", f))
  cat("Get modularity for", tissue, "\n")
  modularity=readRDS(f)
  mod=c(mod, modularity$modularity[length(modularity$modularity)])
  names(mod)[length(mod)]=tissue
  clusters=c(clusters, nrow(modularity$Qcoms))
  names(clusters)[length(clusters)]=tissue
}
saveRDS(mod, file=paste0(opt$directory, "modularity.rds"))
write.table(data.frame("Tissue"=gsub("_", "", names(mod)), "Modularity"=mod, "Clusters"=clusters),
            col.names = TRUE, row.names = FALSE, quote=F, sep="\t", file=paste0(opt$directory, "modularity.txt"))

