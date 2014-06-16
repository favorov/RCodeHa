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

#we output: DM.M.noodles.bonf.bed - for Bonferroni-corrected p-val<0.05
#seqnames start end width strand p.value ishyper

#DM.M.noodles.bonf.strict.bed - pure bed with the DM intervals for GREAT analysis

#DM.M.noodles.bonf.adjacent.genes.tsv genes overlapped (flanks, possibly) by DM noodles
#seqnames	start	end	width	strand	gene_id	SYMBOL	ishyper	p.value

#and, the same 3 files with fdr<0.1 instead of Bonferroni correction: 
#DM.M.noodles.fdr.bed, DM.M.noodles.fdr.strict.bed, DM.M.noodles.fdr.adjacent.genes.tsv

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
	stop('cannot load correct noodles.M.Rda')
	

noodles.M.fisher.results.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('noodles.M.fisher.results.Rda'))
	if ('fisher.noodles.M.result' %in% load('noodles.M.fisher.results.Rda'))
		if (class(fisher.noodles.M.result)=='data.frame')
			if('fisher.p.values' %in% colnames(fisher.noodles.M.result))
				noodles.M.fisher.results.loaded<-TRUE

if(!noodles.M.fisher.results.loaded)
	stop('cannot load correct noodles.M.fisher.results.Rda')

message('all loaded')


#bonferroni
message('bonferroni')
DM.M.noodles.indices<-which(p.adjust(fisher.noodles.M.result$fisher.p.values,method='bonferroni')<0.05)
DM.M.noodles<-noodles.M[DM.M.noodles.indices]
DM.M.noodles$p.value<-fisher.noodles.M.result$fisher.p.values[DM.M.noodles.indices]
DM.M.noodles$ishyper<-fisher.noodles.M.result$CI_95_L[DM.M.noodles.indices]>1

DM.Genes<-gene.list.by.overlap(noodles=DM.M.noodles,flanks=flanks)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.M.noodles.bonf.adjacent.genes.tsv',sep='\t',row.names=FALSE,quote=FALSE)

DM.M.table<-as(DM.M.noodles,'data.frame')

write.table(DM.M.table,file='DM.M.noodles.bonf.bed',quote=FALSE,row.names=FALSE)

con<-file('DM.M.noodles.bonf.strict.bed')
export(DM.M.noodles,con,'bed',ignore.strand=TRUE)
close(con)

#fdr
message('fdr')
DM.M.noodles.indices<-which(p.adjust(fisher.noodles.M.result$fisher.p.values,method='fdr')<0.1)
DM.M.noodles<-noodles.M[DM.M.noodles.indices]
DM.M.noodles$p.value<-fisher.noodles.M.result$fisher.p.values[DM.M.noodles.indices]
DM.M.noodles$ishyper<-fisher.noodles.M.result$CI_95_L[DM.M.noodles.indices]>1

DM.Genes<-gene.list.by.overlap(noodles=DM.M.noodles,flanks=flanks)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.M.noodles.fdr.adjacent.genes.tsv',sep='\t',row.names=FALSE,quote=FALSE)

DM.M.table<-as(DM.M.noodles,'data.frame')

write.table(DM.M.table,file='DM.M.noodles.fdr.bed',quote=FALSE,row.names=FALSE)

con<-file('DM.M.noodles.fdr.strict.bed')
export(DM.M.noodles,con,'bed',ignore.strand=TRUE)
close(con)

