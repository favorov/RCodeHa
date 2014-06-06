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

DM.Genes<-gene.list.by.overlap(noodles=DM.M.noodles,flanks=flanks)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.Genes.short.by.M.noodles.tsv',sep='\t',row.names=FALSE)

