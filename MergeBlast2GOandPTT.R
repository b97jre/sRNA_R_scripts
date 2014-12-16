# Functions for remapping the geneannotation from PID to Spy genes

setwd("/Users/johanreimegard/Vetenskap/Data/sRNA_Charpentier/blast2GO/blast2go")

pptSystem = read.table("NC_002737.ver2.txt",header=TRUE, sep ="\t")
tail(pptSystem)

GOterms = read.table("annot_GOs_20140827_1343.txt" ,sep = "\t", header = FALSE,quote="")

dim(GOterms)

colnames(GOterms) <- c("PID","Name","GOterm","Function","GOdomain") 

MergedInfo  = merge(x=GOterms,y=pptSystem,by="PID")


head(MergedInfo)
NewGOtermName <- MergedInfo[,c(10,9,3,4,5)]

head(NewGOtermName)
write.table(NewGOtermName,file="NC_002737_GOtermAnnotation.tab.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(MergedInfo,file="NC_002737_GOtermAnnotation.tab.extraInfo.txt",quote=FALSE,sep="\t",row.names=FALSE)



