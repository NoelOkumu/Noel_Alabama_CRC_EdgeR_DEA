#This file is to be used for DEA using EdgeR

#Set working Directory
Alabama="/home/noel/Noel_Meru_Server/HTSeq_Counts"
setwd(Alabama)

#Read in Gene mapping
mapping<- read.table("~/Noel_Meru_Server/HTSeq_Counts/ENSG_ID2Name.txt", header= FALSE,stringsAsFactors = FALSE,row.names=1)

#Read in the countsmatrix
Alabama_rawdata_count<- read.table("~/Noel_Meru_Server/HTSeq_Counts/Alabama_HTSeq_Counts.tsv", header=TRUE, stringsAsFactors = FALSE,row.names=1)

#Check dimensions
dim(Alabama_rawdata_count)

#Filtering
#Require atleast 1/12 of samples to have expressed count >=10
sample_cutoff<- (1/12)
count_cutoff<-5

#Define a function to calculate the fraction of values
 getFE<-function(data, count_cutoff){
     FE<-(sum(data>=count_cutoff) / length(data))
    return (FE)
   }

#Apply the function to all genes and filter out genes not meeting the sample cutoff
fraction_expressed<-apply(Alabama_rawdata_count,1,getFE,count_cutoff)
keep<-which(fraction_expressed>=sample_cutoff)
Alabama_rawdata_count<-Alabama_rawdata_count[keep, ]

#Check dimensions to see filtering effect
dim(Alabama_rawdata_count)

################
##Run EdgeR##
################

#Import tool
library("edgeR")
library("edgeR")

#Make class labels to be used in EdgeR
class <- paste(rep(c("UMCV22", "UMCV23", "UMCV24", "UMCV25", "UMCV27", "UMCV30"), each = 2), c("_N", "_T"), sep = "")
class<-as.factor(class)

#Get common gene names
Gene<-rownames(Alabama_rawdata_count)
Symbol<-mapping[Gene,1]
gene_annotations<-cbind(Gene,Symbol)

#Make DGElist Object
y<-DGEList(counts=rawdata,genes=gene_annotations,group=class)

#Normalise the entire dataset
y<-calcNormFactors(y)

#Estimate dispersion on entire dataset
#Apply estimateCommonDisp and estimateTagwiseDisp on y
y<-estimateCommonDisp(y)
y<estimateTagwiseDisp(y)
y
