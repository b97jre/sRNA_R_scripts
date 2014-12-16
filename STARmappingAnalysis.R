
library("ggplot2")
library("plyr")

setwd("/Users/johanreimegard/Vetenskap/Data/sRNA_Charpentier/2014")

mappingInfo <- read.table("Log.summary.ver2",sep="\t")

colNames <- c("Started_job_on","Started_mapping_on","Finished_on","Mapping_speed_,_Million_of_reads_per_hour","Number_of_input_reads","Average_input_read_length","Uniquely_mapped_reads","Uniquely_mapped_reads_Percent","Average_mapped_length","Number_of_splices_Total","Number_of_splices_Annotated_(sjdb)","Number_of_splices_GT/AG","Number_of_splices_GC/AG","Number_of_splices_AT/AC","Number_of_splices_Non-canonical","MismatchRate","Deletion_rate_per_base","Deletion_average_length","Insertion_rate_per_base","Insertion_average_length","Number_of_reads_mapped_to_multiple_loci","Percent_of_reads_mapped_to_multiple_loci","Number_of_reads_mapped_to_too_many_loci","Percent_of_reads_mapped_to_too_many_loci","Percent_of_reads_unmapped_too_many_mismatches","Percent_of_reads_unmapped_too_short","Percent_of_reads_unmapped_other","fileName")

colnames(mappingInfo)<-colNames 

test <-strsplit(as.character(mappingInfo$fileName),split="\\.")
df <- ldply(test)
colnames(df) <- c("seqRound", "Sample_Reference","log")
mappingInfo$seqRound = df$seqRound 

test <-strsplit(as.character(df$Sample_Reference),split="_R_")
df <- ldply(test)
colnames(df) <- c("Sample", "Reference")

mappingInfo$Sample = df$Sample 
mappingInfo$ReadType = "NO_QC_Paired_End" 

mappingInfo$seqRound[mappingInfo$seqRound == "STAR_untreated"] = "FIRST_ROUND"
mappingInfo$seqRound[mappingInfo$seqRound == "STAR_untreated2"] = "SECOND_ROUND"




mappingInfo2 <- read.table("TrimmedReads.STARsummary.txt",sep="\t")

colNames <- c("Started_job_on","Started_mapping_on","Finished_on","Mapping_speed_,_Million_of_reads_per_hour","Number_of_input_reads","Average_input_read_length","Uniquely_mapped_reads","Uniquely_mapped_reads_Percent","Average_mapped_length","Number_of_splices_Total","Number_of_splices_Annotated_(sjdb)","Number_of_splices_GT/AG","Number_of_splices_GC/AG","Number_of_splices_AT/AC","Number_of_splices_Non-canonical","MismatchRate","Deletion_rate_per_base","Deletion_average_length","Insertion_rate_per_base","Insertion_average_length","Number_of_reads_mapped_to_multiple_loci","Percent_of_reads_mapped_to_multiple_loci","Number_of_reads_mapped_to_too_many_loci","Percent_of_reads_mapped_to_too_many_loci","Percent_of_reads_unmapped_too_many_mismatches","Percent_of_reads_unmapped_too_short","Percent_of_reads_unmapped_other","fileName")
colnames(mappingInfo2)<-colNames 


mappingInfoMerged = mappingInfo2[grepl("merged",mappingInfo2$fileName),]

test <-strsplit(as.character(mappingInfoMerged$fileName),split="\\.")
df <- ldply(test)
colnames(df) <- c("seqRound", "Sample","log")
mappingInfoMerged$seqRound = df$seqRound 
mappingInfoMerged$Sample = df$Sample 
mappingInfoMerged$ReadType = "QC_Merged" 


mappingInfoOther = mappingInfo2[!grepl("merged",mappingInfo2$fileName),]
test <-strsplit(as.character(mappingInfoOther$fileName),split="\\.")
df <- ldply(test)
colnames(df) <- c("seqRound", "Sample","log")
mappingInfoOther$seqRound = df$seqRound 
mappingInfoOther$Sample = df$Sample 
mappingInfoOther$ReadType = "QC_PairedEnd" 


mappingInfo = rbind(mappingInfoMerged,mappingInfoOther,mappingInfo)



test <-strsplit(as.character(mappingInfo$Sample),split="_R")
df <- ldply(test)
colnames(df) <- c("Sample")
mappingInfo$Sample = df$Sample 




ER <- strsplit(as.character(mappingInfo$MismatchRate),split="%")
df <- ldply(ER)
colnames(df) <- c("ErrorRate")
mappingInfo$MismatchRate = as.numeric(df$ErrorRate) 




ER <- strsplit(as.character(mappingInfo$Percent_of_reads_unmapped_too_short),split="%")
df <- ldply(ER)
colnames(df) <- c("ErrorRate")
mappingInfo$Percent_of_reads_unmapped_too_short = as.numeric(df$ErrorRate) 



mappingInfo$fractionMappedReads =mappingInfo$Uniquely_mapped_reads /  mappingInfo$Number_of_input_reads
mappingInfo$fractionMappedReadsMultiple =mappingInfo$Number_of_reads_mapped_to_multiple_loci /  mappingInfo$Number_of_input_reads
mappingInfo$fractionMappedReadsMultiple2 =mappingInfo$Number_of_reads_mapped_to_too_many_loci /  mappingInfo$Number_of_input_reads

mappingInfo$totalMapped =(mappingInfo$Number_of_reads_mapped_to_too_many_loci +  mappingInfo$Uniquely_mapped_reads + mappingInfo$Number_of_reads_mapped_to_multiple_loci)/  mappingInfo$Number_of_input_reads


ggplot(mappingInfo, aes(x=factor(Sample),fill=factor(seqRound),weight=MismatchRate))+ 
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Percent error per base")+xlab("Sample")+ facet_grid(. ~ ReadType)
ggsave(filename="errorRate.pdf")


ggplot(mappingInfo, aes(x=factor(Sample),fill=factor(seqRound),weight=Uniquely_mapped_reads))+ 
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Uniquely mapped reads")+xlab("Sample")+ facet_grid(. ~ ReadType)
ggsave(filename="NrOfMappedReads.pdf")


ggplot(mappingInfo, aes(x=factor(Sample),fill=factor(seqRound),weight=fractionMappedReads))+ 
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Fraction Uniquely Mapped Reads")+xlab("Sample")+ facet_grid(. ~ ReadType)
ggsave(filename="FractionUniquleyMappedReads.pdf")

ggplot(mappingInfo, aes(x=factor(Sample),fill=factor(seqRound),weight=fractionMappedReadsMultiple))+ 
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Fraction Mapped Reads to multiple locations (<=3)")+xlab("Sample")+ facet_grid(. ~ ReadType)
ggsave(filename="FractionMappedReadsMultiple.pdf")

ggplot(mappingInfo, aes(x=factor(Sample),fill=factor(seqRound),weight=fractionMappedReadsMultiple2))+ 
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Fraction Mapped Reads to multiple locations (>3)")+xlab("Sample")+ facet_grid(. ~ ReadType)
ggsave(filename="FractionMappedReadsToomany.pdf")

ggplot(mappingInfo, aes(x=factor(Sample),fill=factor(seqRound),weight=totalMapped))+ 
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Reads mapped (total)")+xlab("Sample")+ facet_grid(. ~ ReadType)
ggsave(filename="FractionMappedReads.pdf")




ggplot(mappingInfo, aes(x=factor(Sample),fill=factor(seqRound),weight=Number_of_input_reads))+ 
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Number of input Reads")+xlab("Sample")+ facet_grid(. ~ ReadType)

ggsave(filename="NrOfReads.pdf")




ggplot(mappingInfo, aes(x=factor(Sample),fill=factor(seqRound),weight=Percent_of_reads_unmapped_too_short))+ 
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Percent")+xlab("Sample")+ facet_grid(. ~ ReadType)
ggsave(filename="PercentToShort.pdf")







