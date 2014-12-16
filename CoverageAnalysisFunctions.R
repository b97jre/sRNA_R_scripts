
library(quantmod)

GetWigInfo <- function(fileName, outFileNamePrefix){
  # Read table with all info
  Distribution <- read.table(fileName,sep="\t", header=TRUE)
  # Change colnames to get more understandable colnames
  Distribution = changeColNames(Coverage=Distribution,from=background, to=background2)
  # Merge 
  Merged_Distribution = mergeColumns(Coverage=Distribution,metaInfoList=metaInfoList,
                                     mergeInfo=names(metaInfoListMerged),
                                     Allinfo=names(metaInfoList))
  MetaInfo_merged = getMetaInfo(Distribution=Merged_Distribution ,metaInfoList=metaInfoListMerged)
  printInfoToFigures(outFileNamePrefix, Merged_Distribution, MetaInfo_merged, maxCoverage, cutoff, genomeName) 
  printInfoToWIGFiles(outFileNamePrefix, Merged_Distribution, MetaInfo_merged, maxCoverage, cutoff,  genomeName = genomeName) 
  
}


mergeColumns <- function(Coverage, metaInfoList,mergeInfo=c("background","treatment"),
                         Allinfo=c("background","treatment","round","readTypes")){
  MetaInfo = getMetaInfo(Distribution=Coverage ,metaInfoList=metaInfoList)
  print(MetaInfo)
  merged_Coverage = get_one_distribution_per_sample(Distribution=Coverage,MetaInfo=MetaInfo,
            mergeInfo=mergeInfo,Allinfo= Allinfo) 
  return (merged_Coverage)
}

changeColNames <- function(Coverage,from,to){
  tempNames = colnames(Coverage)
  for(i in 1:length(from)){
    tempNames = gsub(pattern=from[[i]], replacement=to[[i]], x=tempNames) 
  }
  colnames(Coverage) <- tempNames
  return (Coverage)
}




printInfoToFigures <- function(prefix , Distribution, MetaInfo, maxCoverage = 50, cutoff = 0.05, genomeName ="gi|15674250|ref|NC_002737.1|" , pageInfo= "background"){  
  pdfFileName = paste (prefix,"_Cumulative_Distribution.pdf")
  pdf(pdfFileName)
  print_BED_distribution_one_per_sample(Distribution=Distribution ,MetaInfo=MetaInfo,maxCoverage=maxCoverage)
  dev.off()
  

  LineVariables = as.factor(MetaInfo[[pageInfo]])
  LineVariableLevels = levels(LineVariables)
  
  pdfFileName = paste (prefix,"_Cumulative_Distribution_per_",pageInfo,".pdf")
  pdf(pdfFileName)
  par(mfrow=c(2,ceiling(length(LineVariableLevels)/2)))
  for(i in 1:length(LineVariableLevels)){
    minorDistribution = Distribution[,which(LineVariables == LineVariableLevels[[i]])] 
    minorMetaInfo = MetaInfo[which(LineVariables == LineVariableLevels[[i]]),]
    print_BED_distribution_one_per_sample(Distribution=minorDistribution ,MetaInfo=minorMetaInfo,maxCoverage=maxCoverage, title=LineVariableLevels[[i]])
  }
  dev.off()
  



}

printInfoToWIGFiles <- function(fileName, Distribution, MetaInfo, maxCoverage = 50, cutoff = 0.05, genomeName ="gi|15674250|ref|NC_002737.1|" ){  

  peakTable = findPeaksinData(Distribution=Distribution, cutoff=cutoff,maxCoverage=maxCoverage)
  printPeaksToWig(peakTable=peakTable, prefix=fileName,genomeName=genomeName)
  printPeaksToWig(peakTable=Distribution, prefix=fileName,genomeName=genomeName,viewLimits = paste("0",maxCoverage,sep=":"))
}




get_one_distribution_per_sample <- function(Distribution, MetaInfo,mergeInfo=c("background","treatment"),
                                            Allinfo=c("background","treatment","round","readTypes")){
  MetaInfo$NewColNames = MetaInfo[[Allinfo[[1]]]]
  if(length(Allinfo) > 1){
    for(i in 2:length(Allinfo)){
      MetaInfo$NewColNames = paste(MetaInfo$NewColNames,MetaInfo[[Allinfo[i]]], sep = "_")
    }
  }
  MetaInfo$Experiment = MetaInfo[[mergeInfo[[1]]]]
  if(length(mergeInfo) > 1){
    for(i in 2:length(mergeInfo)){
      MetaInfo$Experiment = paste(MetaInfo$Experiment,MetaInfo[[mergeInfo[i]]], sep = "_")
    }
  }
  MetaInfo$Experiment = as.factor(MetaInfo$Experiment)
  ExprerimentsLevels = levels(MetaInfo$Experiment)
  for(i in 1:length(ExprerimentsLevels)){
    Distribution[ExprerimentsLevels[i]] = rowSums(Distribution[ ,grep(ExprerimentsLevels[i], MetaInfo$NewColNames)])
  }
  Distribution_merged = Distribution[ , colnames(Distribution) %in% ExprerimentsLevels]
  return (Distribution_merged)
}



getMetaInfo <-function(Distribution, metaInfoList){
  MetaInfo = as.data.frame(colnames(Distribution))
  colnames(MetaInfo) <- c("FileNames")
  metaInfoListNames = names(metaInfoList)
  
  for(i in 1:length(metaInfoListNames)){
    MetaInfo[metaInfoListNames[i]] = metaInfoList[metaInfoListNames[i]][[1]][[1]]
    if(length(metaInfoList[metaInfoListNames[i]][[1]])>1){
      for(j in 2:length(metaInfoList[metaInfoListNames[i]][[1]])){
        if(sum(grepl(metaInfoList[metaInfoListNames[i]][[1]][[j]], MetaInfo$FileNames))>0 ){
          MetaInfo[grep(metaInfoList[metaInfoListNames[i]][[1]][[j]], MetaInfo$FileNames),metaInfoListNames[i]] =metaInfoList[metaInfoListNames[i]][[1]][[j]] 
        }
      }
    }
  }
  return(MetaInfo)
  
}


print_BED_distribution_one_per_bamFile <- function(Distribution, Metainfo , geneticBackground ="RNAseYmut" ){
  
  variables = as.factor(MetaInfo$Treatment)
  VariableLevels = levels(variables)
  
  Palete = rainbow(nlevels(variables), alpha = 1)
  PaleteLevels = VariableLevels
  
  LineVariables = as.factor(MetaInfo$readType)
  LineVariableLevels = levels(LineVariables)
  
  PointVariables = as.factor(MetaInfo$round)
  PointVariableLevels = levels(PointVariables)
  
  
  m <- as.matrix(Distribution)
  m[m>20] <- 20
  
  first = 1
  for (i in which(MetaInfo$Background == Ba)){
    data = m[,i]
    breaks = seq(0, 20, by=1)
    duration.cut = cut(data, breaks, right=FALSE)
    duration.freq = table(duration.cut)
    cumfreq0 = c(0, cumsum(duration.freq)/sum(duration.freq))
    
    SpecificColor = Palete[match(variables[i],PaleteLevels)]
    linetype = match(LineVariables[i],LineVariableLevels)
    pointType = match(PointVariables[i],PointVariableLevels)
    
    if(first == 1){
      
      plot(breaks, cumfreq0, col=SpecificColor,ylim=c(0.95,1),lty= linetype,pch = pointType)
      info = paste(variables[i],LineVariables[i],PointVariables[i], sep = ",")
      LineTypes = linetype
      PointTypes = pointType
      ColType = SpecificColor
      first = 2
    }else{
      points(breaks, cumfreq0, col=SpecificColor,lty= linetype,pch = pointType )
      info = c(info, paste(variables[i],LineVariables[i],PointVariables[i], sep = ","))
      LineTypes = c(LineTypes,linetype)
      PointTypes = c(PointTypes,pointType)
      ColType = c(ColType,SpecificColor)
    }
    lines(breaks, cumfreq0, col=SpecificColor,lty= linetype)
  }
  
  legend(x="bottom",legend=VariableLevels,col=Palete, lty = 1)
  legend(x="bottomright",legend = LineVariableLevels,col ="black", lty = 1:length(levels(LineVariables)))
  legend(x="right",legend = PointVariableLevels,col ="black", pch = 1:length(levels(PointVariables)))
  
}



print_BED_distribution_one_per_sample <- function(Distribution, MetaInfo, colInfo = "background" , lineInfo = "treatment", 
                                                  maxCoverage = 100, step = 1, 
                                                  title = "Cumulative Distribution over all positions"){
  
  
  #colourInfo 
  
  variables = as.factor(MetaInfo[[colInfo]])
  VariableLevels = levels(variables)
  if(length(VariableLevels) > 1){
    Palete = rainbow(nlevels(variables), alpha = 1)
  }else{
    Palete = c("black")
  }
  PaleteLevels = VariableLevels
    
  #lineInfo 
  LineVariables = as.factor(MetaInfo[[lineInfo]])
  LineVariableLevels = levels(LineVariables)
  
  
  m <- as.matrix(Distribution)
  m[m>maxCoverage] <- maxCoverage-step
  
  first = 1
  for (i in 1:dim(m)[2]){
    data = m[,i]
    breaks = seq(0, maxCoverage, by=step)
    duration.cut = cut(data, breaks, right=FALSE)
    duration.freq = table(duration.cut)
    cumfreq0 = c(0, cumsum(duration.freq)/sum(duration.freq))
    
    SpecificColor = Palete[match(variables[i],PaleteLevels)]
    #linetype = match(LineVariables[i],LineVariableLevels)
    linetype = 1
    pointType = match(LineVariables[i],LineVariableLevels)
    
    if(first == 1){
      plot(breaks, cumfreq0, col=SpecificColor,ylim=c(0.95,1),lty= linetype,pch = pointType,xlab="Nr of reads",ylab="Cumlative frequency")
      info = paste(variables[i],LineVariables[i], sep = ",")
      LineTypes = linetype
      PointTypes = pointType
      ColType = SpecificColor
      first = 2
    }else{
      points(breaks, cumfreq0, col=SpecificColor,lty= linetype,pch = pointType )
      info = c(info, paste(variables[i],LineVariables[i], sep = ","))
      LineTypes = c(LineTypes,linetype)
      PointTypes = c(PointTypes,pointType)
      ColType = c(ColType,SpecificColor)
    }
    lines(breaks, cumfreq0, col=SpecificColor,lty= linetype)
  }
  
  legend(x="bottomright",legend = info,col =ColType, lty = LineTypes,pch=PointTypes,cex=0.6)
  title(title)
}  


findPeaksinData <- function(Distribution, cutoff = 0.05, maxCoverage = 1000){
  
  distMatrix <- as.matrix(Distribution)
  m = distMatrix
  m[m>maxCoverage] = maxCoverage
  peaks = m
  peaks[peaks > 0 ] = 0 
  
  for (i in 1:dim(m)[2]){
    data = m[,i]
    
    breaks = seq(0, maxCoverage, by=1)
    duration.cut = cut(data, breaks, right=FALSE)
    duration.freq = table(duration.cut)
    cumfreq0 = c(0, cumsum(duration.freq)/sum(duration.freq))
    threshold = sum((1-cumfreq0)>cutoff)
    peaks[findPeaks(data,thresh=threshold),i]=1
  }
  peaks2 <- as.data.frame(peaks)
  return (peaks2)
}  



printPeaksToWig <- function(prefix, genomeName = "gi|15674250|ref|NC_002737.1|", peakTable,viewLimits="-1:1"){
  fileNames = colnames(peakTable)
  for(i in 1:length(fileNames)){
    fileInfo =  paste(prefix, fileNames[i],sep=".")
    fileName = paste(fileInfo,"wig",sep=".")
    firstRow = paste("track name=\"",fileInfo,"\" color=255,0,255 altColor=255,0,255 graphType=bar viewLimits=",viewLimits, sep = "")
    write.table(x=firstRow,file=fileName,quote=FALSE,row.names=FALSE,col.names=FALSE)
    secondRow = paste("fixedStep chrom=",genomeName," start=1 step=1", sep = "")  
    write.table(x=secondRow,file=fileName,quote=FALSE,row.names=FALSE,col.names=FALSE, append = TRUE)
    write.table(x=peakTable[[fileNames[i]]],file=fileName,quote=FALSE,row.names=FALSE,col.names=FALSE, append = TRUE)
    
  }
  
}


