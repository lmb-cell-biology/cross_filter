#!/usr/bin/R

args <- commandArgs(trailingOnly=TRUE)

coverage_list <- args[1]
heatmap_file <- paste0(args[2],"_heatmap.pdf")

library(pheatmap)

all_strains_coverage <- read.table(file = coverage_list,sep = "\t",header=TRUE)

head(all_strains_coverage)

# all_strains_coverage <- data.frame(Genome_cov=rnorm(n = 50,mean = 40,2),Exon_cov=rnorm(n = 50,mean = 40,2))
# row.names(all_strains_coverage)<-paste0("AX65",1:50)

pdf(heatmap_file)
pheatmap(all_strains_coverage,cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE)
dev.off()
