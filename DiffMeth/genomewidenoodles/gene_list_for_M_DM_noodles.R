if (!suppressWarnings(require('Differential.Coverage')))
{
	if (!suppressWarnings(require('devtools')))
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("devtools")
		library("devtools")
	}
	install_github('Differential.Coverage','favorov')
	#load_all('../../../../differential.coverage/')
	library('Differential.Coverage')
}

if (!require('Homo.sapiens'))
{
  source("http://bioconductor.org/biocLite.R")
  biocLite("Homo.sapiens")
  library('Homo.sapiens')  
}
if (!require('org.Hs.eg.db'))
{
  source("http://bioconductor.org/biocLite.R")
  biocLite('org.Hs.eg.db')
  library('org.Hs.eg.db')  
}

flanks<-10000

noodles.M.loaded<-FALSE
# we can the whole thing to noodles.M.Rda
if(file.exists('noodles.M.Rda'))
{
	loaded<-load('noodles.M.Rda')
	if ('noodles.M.methylation' %in% loaded) 
		if (class(noodles.M.methylation)=='data.frame')
			if ('noodles.M' %in% loaded)
				if(class(noodles.M)=='GRanges')
			noodles.M.loaded<-TRUE
}
if(!noodles.M.loaded)
{
	source('noodles_M_gw.R')
}

noodles.M.fisher.results.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('noodles.M.fisher.results.Rda'))
	if ('fisher.p.values' %in% load('noodles.M.fisher.results.Rda'))
			noodles.M.fisher.results.loaded<-TRUE

if(!noodles.M.fisher.results.loaded)
{
	source('noodles_M_gw.R')
}

message('all loaded')

DM.M.noodles.indices<-which(fisher.p.values*tests.number<0.05)
DM.M.noodles<-noodles.M[DM.M.noodles.indices]
DM.M.noodles$p.value<-fisher.p.values[DM.M.noodles.indices]
DM.M.noodles$ishyper<-CI_95_L[DM.M.noodles.indices]>1

expanded.DM.M.noodles<-DM.M.noodles
#inflate DM noodles
start(expanded.DM.M.noodles)<-pmax(0,start(DM.M.noodles)-flanks)
end(expanded.DM.M.noodles)<-pmin(end(DM.M.noodles)+flanks,as.integer(seqlengths(DM.M.noodles)[as.character(seqnames(DM.M.noodles))]))

#prepare gene TSS
TSS<- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

geneSymbols <- select(
	org.Hs.eg.db,
	keys=as.character(TSS$gene_id),
	columns=c('SYMBOL'),
	keytype='ENTREZID'
)

TSS$SYMBOL <- geneSymbols$SYMBOL

tss.start<-ifelse(strand(TSS)=='+',start(TSS),end(TSS))

start(TSS)<-tss.start
end(TSS)<-tss.start

#make overlap

