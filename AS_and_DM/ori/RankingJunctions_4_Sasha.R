##### Make a unique function for both analysis. 
### Only use EVA not DEVas. It might outperform the existing combination 

#juncExp: juncExpr
#juncList: list of the junctions
#
library('Homo.sapiens')
library('org.Hs.eg.db')
library('GenomicRanges')
library("GSReg")
library(EBSeq)
library(limma)


source("C:/Users/bahman/Dropbox/SEVApaper/PaperSuppl/Scripts/functions.R")


#loading gene expression
load("C:/Users/bahman/Dropbox/HaJunctions20150707/25Jul2015_FinalProgressReportAnalysis/RawDat/junc.RPM.rda")
#loading junction expression and phenotype
#load("C:/Postdoc Reserach/Splice Variance/From Theresa/JuncRPM_11Nov2014.RDa")
#phenotype 
NormalSamp <- colnames(junc.RPM)[grep(pattern = "Normal",colnames(junc.RPM))]
TumorSamp <- colnames(junc.RPM)[grep(pattern = "HN",colnames(junc.RPM))]
phenoVect <- c(rep(x= "Normal",length(NormalSamp)),rep(x="Tumor",length(TumorSamp)))
names(phenoVect) <- c(NormalSamp,TumorSamp)

#gene exp removing duplicated names



gn <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
gSymbol <- select(org.Hs.eg.db,keys=as.character(gn$gene_id),
                  columns=c('SYMBOL'),keytype='ENTREZID')
gn$SYMBOL <- gSymbol$SYMBOL

z <- strsplit(sub('-',':',rownames(junc.RPM)),':')
mychr <- sapply(X=z,FUN = function(x) x[1])
mystart <- sapply(X=z,FUN = function(x) x[2])
myend <- sapply(X=z,FUN = function(x) x[3])

### puting junction in GRanges format
juncRanges <- as.data.frame(strsplit(sub('-',':',rownames(junc.RPM)),':'),stringsAsFactors = FALSE)
junctionsGRanges <- GRanges(seqnames = Rle(mychr), 
                            ranges = IRanges(start = as.numeric(mystart), end = as.numeric(myend)))

### finding junction for each gene
overlapJunction <- findOverlaps(junctionsGRanges,gn)

#making overlap matrix
juncnames <- rownames(junc.RPM)
genesJunction <- tapply(rownames(junc.RPM)[queryHits(overlapJunction)],
                        gn$SYMBOL[subjectHits(overlapJunction)],function(x) x)

genesJunctionInd <- tapply(queryHits(overlapJunction),
                           gn$SYMBOL[subjectHits(overlapJunction)],function(x) x)



genesDSxls <-  read.csv("C:/Users/bahman/Dropbox/HaJunctions20150707/25Jul2015_FinalProgressReportAnalysis/Results/EVAGenes.csv")
pvaluesEVAOut <- genesDSxls[,2:3]
row.names(pvaluesEVAOut) <- genesDSxls[[1]]
#pvalueMax <- apply(pvaluesEVAOut,MARGIN = 1,FUN = function(x) ifelse(sum(is.na(x))>0,yes = 1,no = max(x)))
#DSgenesEVAOut <- names(which(pvalueMax<0.01))

pvalue <- as.matrix(na.omit(pvaluesEVAOut[1]))
DSgenesEVAOut <- rownames(pvalue)[which(pvalue<0.01)]



MyRest <- sapply(X = genesJunctionInd[DSgenesEVAOut], function(x) { 
  y <- junctionsGRanges[x];
  w <- findOverlaps(y,y);
  mymat <- matrix(0,nrow = length(y), ncol = length(y),dimnames = list(juncnames[x],juncnames[x]));
  mymat[queryHits(w)+ length(x)*(subjectHits(w)-1)] <- 1
  return(mymat)
})

MyRest <- MyRest[names(which(sapply(MyRest, length)>0))]

junctionRanks <- vector(mode = "list",length=length(MyRest))
names(junctionRanks) <- names(MyRest)

for( i in 1:length(MyRest)){
    juncsingenes <- rownames(MyRest[[i]])
    
    #junctionRanks[[i]] <- matrix(data = 0,nrow = length(juncsingenes), 
    #       ncol = ncol(junc.RPM),dimnames = list(juncsingenes,colnames(junc.RPM)))
    
    junctionRanks[[i]] <- t(sapply(juncsingenes, FUN = function(x)
        {
          z <- names(which(MyRest[[i]][x,]>0))
          
          if(is.vector(junc.RPM[z,])){
            return(vector( mode= "numeric",length = ncol(junc.RPM) ) )
          }else{
            #junctionRanks[[i]][j,] <- apply(junc.RPM[z,],MARGIN = 2,FUN = rank)[j,]
            return(apply(junc.RPM[z,],MARGIN = 2,FUN = rank)[x,])
        }
      }))
#     for( j in juncsingenes){
#       z <- names(which(MyRest[[i]][j,]>0))
# #       for( k in ncol(junc.RPM))
# #         junctionRanks[[i]][j,k] <-  0.5*(mean(junc.RPM[j ,k] < junc.RPM[z ,k]) 
# #                                          -  mean(junc.RPM[j ,k] > junc.RPM[z ,k]))+0.5
#       if(is.vector(junc.RPM[z,])){
#         junctionRanks[[i]][j,] <- 0
#       }else{
#         junctionRanks[[i]][j,] <- apply(junc.RPM[z,],MARGIN = 2,FUN = rank)[j,]
#       }
#     }
    if(i%%20 == 0 )
      print(i)

}

save(list = "junctionRanks", file = 
      "C:/Users/bahman/Dropbox/HaJunctions20150707/ForSasha/RankingJunctions_Sasha.rda")
