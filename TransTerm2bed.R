setwd("/Users/johanreimegard/Vetenskap/Data/sRNA_Charpentier/RhoIndependentTerminators")

RIT = read.table("NC_002737.transterm.filtered.tab")

genome = "gi|15674250|ref|NC_002737.1|"
color = "207,3,124"

# removing uninteresting column
RIT = RIT[,c(1,2,4,5,6,7,8,9)]
RIT2 = RIT

colnames(RIT) <- c("Name","start","stop","dir","loc","conf","hp", "tail")

RITBED = convertToBED(RIT,genome, color)
RITBEDHighConf = RITBED[RITBED$conf > 75, ]

BEDinfo = RITBEDHighConf
writeBEDToFile(RITBEDHighConf,BEDfile = "NC_002737.transterm.bed")


writeBEDToFile <- function(BEDinfo, BEDfile="BEDfile.bed", trackName = "RhoIndependent Terminators",trackInfo="Predicted by transterm_hp_v2.09" ,itemRGB = "on"){


  FirstString =paste("track name=",dQuote(trackName)," description=",dQuote(trackInfo)," visibility=2 useScore=1" , sep = "")
  write.table(x = FirstString,file =  BEDfile , quote = FALSE,row.names = FALSE,col.names = FALSE)
  write.table(x = BEDinfo,file =  BEDfile , quote = FALSE,row.names = FALSE,col.names = FALSE,append =TRUE, sep = "\t")
  
  
  
}



convertToBED <- function(RIT2, genome ="gi|15674250|ref|NC_002737.1|",color = "207,3,124"){
  RIT2$genome = genome
  RIT2$color = color
  RIT2$left =RIT2$start 
  RIT2$right =RIT2$stop 
  RIT2$left[which(RIT2$stop < RIT2$start)] = RIT2$stop[RIT2$stop < RIT2$start] 
  RIT2$right[which(RIT2$stop < RIT2$start)] = RIT2$start[RIT2$stop < RIT2$start] 
  RIT2$conf = RIT2$conf*10-1 
  
  RITBED <- RIT2[,c("genome","left","right","Name","conf","dir") ]
  head(n = 100, RITBED)
  return (RITBED)
}






colnames(RIT) <- 
head(RIT)


