---
title: "Create LOLA peakome"
output: html_document
---

```{r setup, include=FALSE}
library(GenomicRanges)
library(data.table)
library(dplyr)
library(ggplot2)
library(gridGraphics)
library(splitstackshape)

knitr::opts_chunk$set(echo = TRUE)

codex_dir="/Users/javrodher/Work/biodata/LOLA/nm/t1/resources/regions/LOLACore/hg38/codex/regions/"
encode_dir="/Users/javrodher/Work/biodata/LOLA/nm/t1/resources/regions/LOLACore/hg38/encode_tfbs/regions/"
outdir="/Users/javrodher/Work/RStudio-PRJs/Fatemeh_Hernandolab//results/"
```

```{r}
f_codex=list.files(codex_dir,pattern = "CTCF")
f_encode=list.files(encode_dir,pattern = "Ctcf")
index_codex_file="/Users/javrodher/Work/biodata/LOLA/nm/t1/resources/regions/LOLACore/hg38/codex/index.txt"
index_encode_file="/Users/javrodher/Work/biodata/LOLA/nm/t1/resources/regions/LOLACore/hg38/encode_tfbs/index.txt"

index_codex=read.delim(index_codex_file)
index_encode=read.delim(index_encode_file)

index_codex_ctcf = index_codex[index_codex$antibody=="CTCF",]
index_encode_ctcf = index_encode[index_encode$antibody=="CTCF",]

index_codex_ctcf$mappingGenome
index_encode_ctcf=index_encode_ctcf[index_encode_ctcf$treatment=="None",]
table(index_encode_ctcf$cellType)

fnames=index_encode_ctcf$filename

X=data.frame(V1=NA,V2=NA,V3=NA,V4=NA,stringsAsFactors = F)
for (fname in fnames){
  #fname=fnames[1]
  
  print(fname)
  x = read.delim(paste0(encode_dir,"/",fname),header = F)
  x$V4=fname
  X=rbind(X,x)
}

X=X[-1,]
table(X$V4)

write.table(X[,1:3],paste0(outdir,"/","ctcf_encode_LOLA_hg38.bed"),quote = F,sep="\t",col.names = F,row.names = F)
```

```{r}
## RUN in terminal
# sort -k1,1 -k2,2n ctcf_encode_LOLA_hg38.bed > ctcf_encode_sorted.bed
# bedtools merge -i ctcf_encode_sorted.bed > ctcf_encode_merged.bed
# wc -l ctcf_encode_merged.bed
```
