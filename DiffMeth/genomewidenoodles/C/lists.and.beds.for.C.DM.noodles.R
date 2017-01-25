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

#we output: DM.C.noodles.bonf.bed - for Bonferroni-corrected p-val<0.05
#seqnames start end width strand p.value ishyper

#DM.C.noodles.bonf.strict.bed - pure bed with the DM intervals for GREAT analysis

#DM.C.noodles.bonf.adjacent.genes.tsv genes overlapped (flanks, possibly) by DM noodles
#seqnames	start	end	width	strand	gene_id	SYMBOL	ishyper	p.value

#and, the same 3 files with fdr<0.1 instead of Bonferroni correction: 
#DM.C.noodles.fdr.bed, DM.C.noodles.fdr.strict.bed, DM.C.noodles.fdr.adjacent.genes.tsv

flanks<-10000

noodles.C.loaded<-FALSE
# we can the whole thing to noodles.C.Rda
if(file.exists('noodles.C.Rda'))
{
	loaded<-load('noodles.C.Rda')
	if ('noodles.C.methylation' %in% loaded) 
		if (class(noodles.C.methylation)=='dgCMatrix')
			if ('noodles.C' %in% loaded)
				if(class(noodles.C)=='GRanges')
			noodles.C.loaded<-TRUE
}
if(!noodles.C.loaded)
	stop('cannot load correct noodles.C.Rda')
	

noodles.C.fisher.results.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('noodles.C.fisher.results.Rda'))
	if ('fisher.noodles.C.result' %in% load('noodles.C.fisher.results.Rda'))
		if (class(fisher.noodles.C.result)=='data.frame')
			if('fisher.p.values' %in% colnames(fisher.noodles.C.result))
				noodles.C.fisher.results.loaded<-TRUE

if(!noodles.C.fisher.results.loaded)
	stop('cannot load correct noodles.C.fisher.results.Rda')

message('all loaded')


#bonferroni
message('bonferroni')
DM.C.noodles.indices<-which(p.adjust(fisher.noodles.C.result$fisher.p.values,method='bonferroni')<0.05)
DM.C.noodles<-noodles.C[DM.C.noodles.indices]
DM.C.noodles$p.value<-fisher.noodles.C.result$fisher.p.values[DM.C.noodles.indices]
DM.C.noodles$ishyper<-fisher.noodles.C.result$CI_95_L[DM.C.noodles.indices]>1

DM.Genes<-genes.with.TSS.covered(noodles=DM.C.noodles,flanks=flanks)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.C.noodles.bonf.adjacent.genes.tsv',sep='\t',row.names=FALSE,quote=FALSE)

DM.C.table<-as(DM.C.noodles,'data.frame')

write.table(DM.C.table,file='DM.C.noodles.bonf.bed',quote=FALSE,row.names=FALSE)

con<-file('DM.C.noodles.bonf.strict.bed')
export(DM.C.noodles,con,'bed',ignore.strand=TRUE)
close(con)

#bonferroni with flanks 100000
message('bonferroni,flanks 100000')

DM.Genes<-genes.with.TSS.covered(noodles=DM.C.noodles,flanks=100000)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.C.noodles.bonf.adjacent.genes.100000.tsv',sep='\t',row.names=FALSE,quote=FALSE)

#bonferroni with flanks 1000000
message('bonferroni,flanks 1000000')

DM.Genes<-genes.with.TSS.covered(noodles=DM.C.noodles,flanks=1000000)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.C.noodles.bonf.adjacent.genes.1000000.tsv',sep='\t',row.names=FALSE,quote=FALSE)

#fdr
message('fdr')
DM.C.noodles.indices<-which(p.adjust(fisher.noodles.C.result$fisher.p.values,method='fdr')<0.1)
DM.C.noodles<-noodles.C[DM.C.noodles.indices]
DM.C.noodles$p.value<-fisher.noodles.C.result$fisher.p.values[DM.C.noodles.indices]
DM.C.noodles$ishyper<-fisher.noodles.C.result$CI_95_L[DM.C.noodles.indices]>1

DM.Genes<-genes.with.TSS.covered(noodles=DM.C.noodles,flanks=flanks)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.C.noodles.fdr.adjacent.genes.tsv',sep='\t',row.names=FALSE,quote=FALSE)

DM.C.table<-as(DM.C.noodles,'data.frame')

write.table(DM.C.table,file='DM.C.noodles.fdr.bed',quote=FALSE,row.names=FALSE)

con<-file('DM.C.noodles.fdr.strict.bed')
export(DM.C.noodles,con,'bed',ignore.strand=TRUE)
close(con)

#fdr with no flanks
message('fdr/no flanks')

DM.Genes<-genes.with.TSS.covered(noodles=DM.C.noodles,flanks=0)

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.C.noodles.fdr.adjacent.genes.noflanks.tsv',sep='\t',row.names=FALSE,quote=FALSE)

#fdr with 100 000 flanks
message('fdr/100 000 flanks')

DM.Genes<-genes.with.TSS.covered(noodles=DM.C.noodles,flanks=100000)

save(file='DM.Genes.FDR.100000.Rda',list=c(DM.genes))

DM.Genes.df<-as(DM.Genes,'data.frame')

DM.Genes.df<-DM.Genes.df[order(DM.Genes.df$seqnames,DM.Genes.df$start),]

write.table(DM.Genes.df,file='DM.C.noodles.fdr.adjacent.genes.100000.tsv',sep='\t',row.names=FALSE,quote=FALSE)


