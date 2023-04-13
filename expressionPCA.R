#################################
# SHELL/R script to generate a PCA plot for gene expression datasets
# David Fournier october 2017
# JGU Mainz, Germany
#################################

##SHELL SCRIPT

# downloading of dataset from example
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47029/matrix/GSE47029_series_matrix.txt.gz
gzip -d GSE47029_series_matrix.txt.gz # unzipping the data

cat GSE47029_series_matrix.txt | grep -v \! > data.txt # extracts the data
cat GSE47029_series_matrix.txt | grep \!Sample_characteristics > mdata.txt # extracts the metadata starting by "!Sample_characteristics" in the text file

## R SCRIPT

library("DESeq")#bioconductor package to perform the differential expression analysis
library("pcaMethods")#bioconductor package to perform PCA

data=read.table("data.txt",h=T,fill=T)
mdata=read.table("mdata.txt",h=F,fill=T)
data=data[,2:ncol(data)]
mdata=mdata[,2:ncol(mdata)]

data=data[complete.cases(data), ] #remove rows with NA values
data[data<0]=0

factors=factor(unlist(mdata[2,])) # selects the desired metadata. values used in the article: 2:strain; 5: Sex; 6: genotype; etc. 
counts=newCountDataSet(round(data), factors)

sizef=estimateSizeFactors(counts)

disps=estimateDispersions(sizef) 
#if no convergence, use alternative method: 
disps=estimateDispersions(sizef ,method="pooled-CR", fitType="local") 

plotDispEsts(disps);

vsdFull=varianceStabilizingTransformation(disps)
print(plotPCA(vsdFull))



