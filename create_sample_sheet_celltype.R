ss_file="sample_sheet_predict.txt"
ss=read.delim(ss_file,header=F)
samples=ss[,2]
samples=gsub(samples,pattern=suffix,replacement = "")
ss=data.frame(sample_name=samples,celltype=NA)
write.csv(ss,"sample_sheet_celltype.csv",row.names = F)
