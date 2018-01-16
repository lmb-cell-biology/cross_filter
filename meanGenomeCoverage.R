#!/usr/bin/R

args <- commandArgs(trailingOnly=TRUE)

if(! "data.table" %in% rownames(installed.packages())){
  cat("data.table has not been installed....\nInstalling data.table\n")
  install.packages("data.table")
}

library(data.table)

genomecov <- fread(args[1],header=FALSE,stringsAsFactors=FALSE)
setnames(genomecov,old=names(genomecov),new=c("chr","depth","bases","ttl_length","fraction"))
mean_genomecov <- sum(genomecov[chr=="genome"][,depth*fraction])
cat(paste0("Genome overall coverage = ",mean_genomecov,"\n"))

exoncov <- fread(args[2],header=FALSE,stringsAsFactors=FALSE)
setnames(exoncov,old=names(exoncov),new=c("chr","depth","bases","ttl_length","fraction"))
mean_exoncov <- sum(exoncov[,depth*fraction])
cat(paste0("Overall coverage of exons in genome = ",mean_exoncov,"\n"))

