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
if ( is.null(opt$output    ) ) {opt$output    = './Results/Figures/'     }

cat("Running compute_stats.R with options:\n",
    "directory =", opt$directory, "\n",
    "output =", opt$output, "\n")

### Create output folder if not exists
if (!(dir.exists(opt$output))){dir.create(opt$output, recursive = T, showWarnings = FALSE)}

### Load data and plot matrix network
files=list.files(path=opt$directory, pattern = '.*.rds')
for (f in paste0(opt$directory, files)){
    tissue=gsub(".*/Edges_matrixeqtl_", "", gsub("_clusters.rds", "", f))
    cat("Get modularity for", tissue, "\n")
    modularity=readRDS(paste0(f))
    modularity$edges$red <- as.factor(modularity$edges$red)
    modularity$edges$blue <- as.factor(modularity$edges$blue)
    modularity$red.memb$red.names <- as.factor(modularity$red.memb$red.names)
    modularity$blue.memb$blue.names <- as.factor(modularity$blue.memb$blue.names)
    
    ## Plot networks as matrix
    plotfile = paste0(opt$output, "Networks_", tissue, ".tiff")
    cols = rep("dodgerblue",length(unique(modularity$red.memb$com)))
    
    tiff(plotfile)
    condor.plot.communities(modularity, color_list=cols)
    dev.off()
}

function (condor.object, color_list, point.size = 0.01, xlab = "SNP", 
          ylab = "Gene") 
{
    condor.object=modularity
    color_list=cols
    point.size = 0.01
    xlab = "SNP"
    ylab = "Gene"
    dt0 <- data.table(condor.object$edges)
    setnames(dt0, 1:2, c("SNP", "gene"))
    dt1 <- data.table(condor.object$red.memb)
    setnames(dt1, c("SNP", "red.memb"))
    dt2 <- data.table(condor.object$blue.memb)
    setnames(dt2, c("gene", "blue.memb"))
    dt3 <- merge(dt0, dt1, by = "SNP", all.x = TRUE)
    eqtl_object <- merge(dt3, dt2, by = "gene", all.x = TRUE)
    setkey(eqtl_object, "SNP")
    eqtl_all <- data.table(eqtl_object[!is.na(SNP)])
    eqtl_block <- eqtl_all[blue.memb == red.memb]
    if (nlevels(as.factor(eqtl_block$SNP)) != length(unique(eqtl_block$SNP))) {
        print("warning: empty levels in SNP column. This may cause silent issues with plotting.")
    }
    setkey(eqtl_block, "blue.memb", "red.memb")
    red_tmp <- data.table(rindx = 1:nlevels(as.factor(eqtl_block$SNP)), 
                          SNP = unique(eqtl_block$SNP))
    red_indx <- merge(red_tmp, unique(eqtl_block, by = "SNP")[, 
                                                              c("SNP", "red.memb"), with = FALSE], by = "SNP")
    red_indx[, `:=`(red.com.size, length(unique(SNP))), by = red.memb]
    red_indx[red.com.size > 1, `:=`(rindx, sample(x = rindx)), 
             by = red.memb][, `:=`(red.memb, NULL), ]
    setkey(red_indx, "SNP")
    blue_tmp <- data.table(bindx = 1:nlevels(as.factor(eqtl_block$gene)), 
                           gene = unique(eqtl_block$gene))
    blue_indx <- merge(blue_tmp, unique(eqtl_block, by = "gene")[, 
                                                                 c("gene", "blue.memb"), with = FALSE], by = "gene")
    blue_indx[, `:=`(blue.com.size, length(unique(gene))), by = blue.memb]
    blue_indx[blue.com.size > 1, `:=`(bindx, sample(x = bindx)), 
              by = blue.memb][, `:=`(blue.memb, NULL), ]
    setkey(blue_indx, "gene")
    if (dim(red_indx)[1] != nlevels(as.factor(eqtl_all$SNP)) && dim(blue_indx)[1] != 
        nlevels(as.factor(eqtl_all$gene))) {
        print("Warning! not all nodes in block!")
    }
    m1 <- merge(eqtl_all, red_indx, by = "SNP", all = TRUE)
    m2 <- merge(m1, blue_indx, by = "gene", all = TRUE)
    par(mar = c(3, 3, 3, 0.5) + 0.1)
    m2[red.memb != blue.memb][plot(rindx, bindx, cex = point.size, 
                                   xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", ylim = c(0, 
                                                                                            max(m2$bindx) + 1), xlim = c(0, max(m2$rindx) + 1), 
                                   xlab = "", ylab = "", pch = 19)]
    m2[red.memb == blue.memb][points(rindx, bindx, cex = point.size, 
                                     pch = 19, col = color_list[red.memb])]
    box(lwd = 2)
    mtext(xlab, side = 3, font = 2, cex = 2.5, padj = -0.25)
    mtext(ylab, side = 2, font = 2, cex = 2.5, padj = -0.5)
    cs <- cumsum(rle(sort(m2[!duplicated(SNP)]$red.memb))$lengths)
    lens <- rle(sort(m2[!duplicated(SNP)]$red.memb))$lengths
    lpts <- cs - lens/2
    axis(1, at = lpts, labels = 1:length(color_list), lwd.ticks = -0.1, 
         cex.axis = 1.25, padj = 0.25, font = 2)
}
