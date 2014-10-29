if (!suppressWarnings(require('differential.coverage')))
{
	if (!suppressWarnings(require('devtools')))
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("devtools")
		library("devtools")
	}
	install_github('favorov/differential.coverage')
	#load_all('../../../../../differential.coverage/')
	library('differential.coverage')
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

DM.Genes<-genes.with.TSS.covered(noodles=DM.M.noodles,flanks=flanks)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.M.noodles.bonf.adjacent.genes.tsv',sep='\t',row.names=FALSE,quote=FALSE)

DM.M.table<-as(DM.M.noodles,'data.frame')

write.table(DM.M.table,file='DM.M.noodles.bonf.bed',quote=FALSE,row.names=FALSE)

con<-file('DM.M.noodles.bonf.strict.bed')
export(DM.M.noodles,con,'bed',ignore.strand=TRUE)
close(con)

#bonferroni with flanks 100000
message('bonferroni,flanks 100000')

DM.Genes<-genes.with.TSS.covered(noodles=DM.M.noodles,flanks=100000)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.M.noodles.bonf.adjacent.genes.100000.tsv',sep='\t',row.names=FALSE,quote=FALSE)

#bonferroni with flanks 1000000
message('bonferroni,flanks 1000000')

DM.Genes<-genes.with.TSS.covered(noodles=DM.M.noodles,flanks=1000000)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.M.noodles.bonf.adjacent.genes.1000000.tsv',sep='\t',row.names=FALSE,quote=FALSE)

#fdr
message('fdr')
DM.M.noodles.indices<-which(p.adjust(fisher.noodles.M.result$fisher.p.values,method='fdr')<0.1)
DM.M.noodles<-noodles.M[DM.M.noodles.indices]
DM.M.noodles$p.value<-fisher.noodles.M.result$fisher.p.values[DM.M.noodles.indices]
DM.M.noodles$ishyper<-fisher.noodles.M.result$CI_95_L[DM.M.noodles.indices]>1

DM.Genes<-genes.with.TSS.covered(noodles=DM.M.noodles,flanks=flanks)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.M.noodles.fdr.adjacent.genes.tsv',sep='\t',row.names=FALSE,quote=FALSE)

DM.M.table<-as(DM.M.noodles,'data.frame')

write.table(DM.M.table,file='DM.M.noodles.fdr.bed',quote=FALSE,row.names=FALSE)

con<-file('DM.M.noodles.fdr.strict.bed')
export(DM.M.noodles,con,'bed',ignore.strand=TRUE)
close(con)

#fdr with no flanks
message('fdr/no flanks')

DM.Genes<-genes.with.TSS.covered(noodles=DM.M.noodles,flanks=0)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.M.noodles.fdr.adjacent.genes.noflanks.tsv',sep='\t',row.names=FALSE,quote=FALSE)

#fdr with 100 000 flanks
message('fdr/100 000 flanks')

DM.Genes<-genes.with.TSS.covered(noodles=DM.M.noodles,flanks=100000)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.M.noodles.fdr.adjacent.genes.100000.tsv',sep='\t',row.names=FALSE,quote=FALSE)


