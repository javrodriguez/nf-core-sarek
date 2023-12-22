# Functions
compilePeaks = function(sample,inpdir){
  x=read.delim(paste0(inpdir,"/",sample,suffix,"/maxatac_predict_32bp.bed"),header = F)
  x=x[order(x$V4,decreasing = T),]
  x$V5=sample
  x$V6=1:nrow(x)
  return(x)
}

getStats=function(predicted_ctcf,ctcf_ref,top_list,ss){
  ctcf_ref.gr=makeGRangesFromDataFrame(ctcf_ref,keep.extra.columns = T)
  DF=data.frame(sample=NA,tps=NA,fps=NA,n_peaks=NA,frac_tps=NA,top=NA,stringsAsFactors = F)
  
  for (top in top_list){
    if(top != "all") { 
      predicted=predicted_ctcf[predicted_ctcf$rank %in% 1:top,] 
      } else {predicted = predicted_ctcf }

    predicted.gr=makeGRangesFromDataFrame(predicted,keep.extra.columns = T)
    ovl=as.data.frame(findOverlaps(predicted.gr,ctcf_ref.gr))
    predicted$overlap_lola=F
    predicted$overlap_lola[unique(ovl$queryHits)]=T
    tps = as.data.frame(table(predicted$sample[predicted$overlap_lola==T]))
    fps = as.data.frame(table(predicted$sample[predicted$overlap_lola==F]))
    df=data.frame(sample=tps$Var1,tps=tps$Freq,fps=fps$Freq,stringsAsFactors = F)
    df$n_peaks=df$tps+df$fps
    df$frac_tps=df$tps/df$n_peaks
    df$top=top
    DF=rbind(DF,df)
}  
  DF=DF[-1,]
  DF$celltype=unlist(lapply(DF$sample,function(x) ss$celltype[ss$sample_name==x])) 
  DF$top=as.factor(DF$top)
  DF$top=factor(DF$top,levels=top_list)
  return(DF)
}

# Load libraries
library(data.table)
library(ggplot2)
library(dplyr)

# Set parameters
inpdir="/Users/javrodher/Work/RStudio-PRJs/Fatemeh_Hernandolab//data/maxATAC/"
ctcf_peaks_file="/Users/javrodher/Work/RStudio-PRJs/Fatemeh_Hernandolab/results/ctcf_encode_merged.bed"
ctcf_motif_file="/Users/javrodher/Work/biodata/genomes/hg38/ctcf_motif_hg38.tsv"
top_list=c("all",seq(20000,5000,-5000))
suffix="-maxATAC-predict"

# Load data
ctcf_lola = read.delim(ctcf_peaks_file,header = F)
ctcf_fimo = fread(ctcf_motif_file,header = T)
ctcf_fimo = as.data.frame(ctcf_fimo[-1,3:5])
names(ctcf_lola)=c("chr","start","end")
names(ctcf_fimo)=c("chr","start","end")
ss=read.csv("sample_sheet_celltype.csv", header = T)

# Compile data
res = lapply(ss$sample_name,compilePeaks,inpdir=inpdir)
predicted_ctcf=do.call(rbind.data.frame,res)
names(predicted_ctcf)=c("chr","start","end","score","sample","rank")
predicted_ctcf$celltype=unlist(lapply(predicted_ctcf$sample,function(x) ss$celltype[ss$sample_name==x])) 
fwrite(predicted_ctcf,"predicted_ctcf_data.csv",row.names=F)

# Get stats for QC
stats_lola=getStats(predicted_ctcf, ctcf_ref=ctcf_lola, top_list, ss)
stats_fimo=getStats(predicted_ctcf, ctcf_ref=ctcf_fimo, top_list, ss)
write.csv(stats_lola,"predicted_ctcf_stats_lola.csv",row.names=F)
write.csv(stats_fimo,"predicted_ctcf_stats_fimo.csv",row.names=F)

# Peak count per sample

pdf(peak_count.pdf)
ggplot(stats_lola[stats_lola$top=="all",],aes(x=sample,y=n_peaks))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))+
  ylab("Predicted CTCF peaks (count)")+
  xlab("")
dev.off()
                                  
# Distribution of peak scores
pdf(peak_scores.pdf)
ggplot(predicted_ctcf,aes(x=score))+
  geom_histogram()+
  xlab("peak score")+
  facet_wrap(~sample)
dev.off()

# True positive fraction (all peaks vs top peaks) - Reference: LOLA peaks
pdf(true_positive_lola.pdf)
ggplot(stats_lola,aes(x=top,y=frac_tps))+
  geom_boxplot()+
  facet_wrap(~celltype)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  ylab("True positives (fraction)")+
  xlab("Peak number (top)")
dev.off()

# True positive fraction (all peaks vs top peaks) - Reference: FIMO motifs
pdf(true_positive_fimo.pdf)
ggplot(stats_fimo,aes(x=top,y=frac_tps))+
  geom_boxplot()+
  facet_wrap(~celltype)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  ylab("True positive (fraction)")+
  xlab("Peak count (top)")
dev.off()
