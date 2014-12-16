
setwd("/Users/johanreimegard/Vetenskap/Data/sRNA_Charpentier/2014/coverage")

source('/Users/johanreimegard/git/sRNAproject/CoverageAnalysisFunctions.R')


# Info about what kind of metainfo that can be extracted from the colnames 
# First is default if none of the ohters can be found in the colnames
# In this case there is an extra step for converting EC number to more understandable backgrounds
background = c( "wt","EC1636","EC1788","EC2058")
background2 = c( "wt","RNAseIIImut","Cas9mut","RNAseYmut")
readTypes = c("pairedEnd","merged")
round = c("FIRST","SECOND")
treatment = c("None","TAP","TAP_TEX")
metaInfoList = list(readTypes,round,treatment,background2)
names(metaInfoList) <- c("readTypes","round","treatment","background")


# Info for which of the metainfo that shoud be kept after they have been merged.
# In this case readTypes and round will be merged
metaInfoListMerged = list(treatment,background2)
names(metaInfoListMerged) <- c("treatment","background")


#Info for output files
cutoff=0.03
maxCoverage=100
genomeName ="gi|15674250|ref|NC_002737.1|"


GetWigInfo("Tend.forward.distribution","Tend.Forward")
GetWigInfo("Tend.reverse.distribution","Tend.Reverse")
GetWigInfo("Fend.forward.distribution","Fend.Forward")
GetWigInfo("Fend.reverse.distribution","Fend.Reverse")
GetWigInfo("Total.forward.distribution","Full.Forward")
GetWigInfo("Total.reverse.distribution","Full.Reverse")


fileName = "Fend.reverse.distribution" 


Distribution <- read.table(fileName,sep="\t", header=TRUE)
# Change colnames to get more understandable colnames
Distribution = changeColNames(Coverage=Distribution,from=background, to=background2)
# Merge 

Merged_Distribution = mergeColumns(Coverage=Distribution,metaInfoList=metaInfoList,
                                   mergeInfo=names(metaInfoListMerged),
                                   Allinfo=names(metaInfoList))
MetaInfo_merged = getMetaInfo(Distribution=Merged_Distribution ,metaInfoList=metaInfoListMerged)

pageInfo = "background"

LineVariables = as.factor(MetaInfo_merged[[pageInfo]])
LineVariableLevels = levels(LineVariables)

i = 2

minorDistribution = Merged_Distribution[,which(LineVariables == LineVariableLevels[[i]])] 
minorMetaInfo = MetaInfo_merged[which(LineVariables == LineVariableLevels[[i]]),]

sumInfo = minorDistribution$None_RNAseIIImut*minorDistribution$TAP_RNAseIIImut
minorDistributionTrimmed = minorDistribution[which(sumInfo > 400 & sumInfo < 100000),]
location = which(sumInfo > 400  & sumInfo < 100000)
logDis = log(minorDistributionTrimmed$None_RNAseIIImut/minorDistributionTrimmed$TAP_RNAseIIImut)

hist(logDis,breaks=100)

logDisTrimmed = logDis[logDis>-2 & logDis<1]
mean = mean(logDis)
sd = sqrt(var(logDis))

pvalues = pnorm(logDis,mean=mean,sd=sd)
sigLocations = location[which (pvalues > 0.999)]

siglogDisTrimmed = 
  sigLocations
test = (logDis - mean)/sd


wilcox.test(test  )



plot(minorDistributionTrimmed$None_RNAseIIImut,minorDistributionTrimmed$TAP_RNAseIIImut)
minorDistribution 

