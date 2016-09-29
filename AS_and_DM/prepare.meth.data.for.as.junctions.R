if (!suppressWarnings(require('differential.coverage')))
{
	if (!suppressWarnings(require('devtools')))
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("devtools")
		library("devtools")
	}
	install_github('favorov/Differential.Coverage')
	#load_all('../../../../differential.coverage/')
	library('differential.coverage')
}

as.and.dm.tables.loaded<-FALSE
if ('mbd.as.junction.coverage' %in% ls() &&
'as.junction.ranges' %in% ls() &&
class(as.junction.ranges)=='GRanges' &&	
(
	class(mbd.as.junction.coverage)=='dgCMatrix' || 
	class(mbd.as.junction.coverage)=='matrix')
)	as.and.dm.tables.loaded <-TRUE
#if they are already in the space, do nothing

if((! as.and.dm.tables.loaded) && file.exists('as.and.dm.tables.Rda'))
{
	message('Loading...')
	loaded<-load('as.and.dm.tables.Rda')
	if (
		'as.junction.ranges' %in% loaded &&
		'mbd.as.junction.coverage' %in%  loaded &&
		class(as.junction.ranges)=='GRanges' &&	
		(
			class(mbd.as.junction.coverage)=='dgCMatrix' || 
			class(mbd.as.junction.coverage)=='matrix')
		)	as.and.dm.tables.loaded <-TRUE
}
#if we loaded it, we do nothing

if(!as.and.dm.tables.loaded)
{
	mbd.bed.dir<-'../../Methylation/bedfiles/'
	chrs<-nucl.chromosomes.hg19()
	mbd.bed.files<-dir(mbd.bed.dir) 
	mbd.bed.files<-mbd.bed.files[grep('All_',mbd.bed.files,invert=TRUE)] # remove two 'All_' files
	mbd.bed.ids<-sapply(strsplit(mbd.bed.files,split='_'),function(x){if(x[2]!='DNA') x[2] else x[3]}) #somhere id in pos 2, somewhere in 3
	mbd.bed.ids<-sub('DNA','',mbd.bed.ids)
	mbd.bed.ids<-sub('PT2','PT',mbd.bed.ids)

	load('RankingJunctions_Sasha.rda')
	splice.bed.ids<-colnames(junctionRanks[[1]])
	splice.bed.ids<-sub('_hg19','',splice.bed.ids)
	splice.bed.ids<-sub('Normal$','',splice.bed.ids)
	splice.bed.ids.norm<-grep('HN',splice.bed.ids,invert = TRUE)
	splice.bed.ids[splice.bed.ids.norm]<-paste0('Normal',splice.bed.ids[splice.bed.ids.norm])
	splice.bed.ids<-sub('^[0-9][0-9]','',splice.bed.ids)

	list.by.gene<-GRangesList(lapply(junctionRanks,function(junction.gene) {
	colnames(junction.gene)<-splice.bed.ids
	juncRanges.gene <- as.data.frame(strsplit(sub('-',':',rownames(junction.gene)),':'),stringsAsFactors = FALSE)
	juncRanges.gene <- data.frame(t(juncRanges.gene),stringsAsFactors=FALSE)
	rownames(juncRanges.gene)<-rownames(junction.gene)
	colnames(juncRanges.gene)<-c('chr','start','end')
	juncRanges.gene<-cbind(juncRanges.gene,junction.gene)
	juncRanges.gene$start<-as.numeric(juncRanges.gene$start)
	juncRanges.gene$end<-as.numeric(juncRanges.gene$end)
	junctionsGRanges.gene <- GRanges(juncRanges.gene,seqinfo=nucl.chromosomes.hg19())
	} ))
	#so, we prepared a juction GRanges object... Now, differential.coverage()
	mcols(list.by.gene)<-data.frame(gene=names(list.by.gene))
	as.junction.ranges<-unlist(list.by.gene)
	mbd.as.junction.coverage<-count.coverage.of.noodles(as.junction.ranges,paste(mbd.bed.dir,mbd.bed.files,sep='/'),mbd.bed.ids)
	as.genes.symbols<-names(junctionRanks)
	save(file='as.and.dm.tables.Rda',list=c('as.junction.ranges','mbd.as.junction.coverage','as.genes.symbols'))
}

noodles.as.genes.tables.loaded<-FALSE
if ('mbd.as.genes.noodles.coverage' %in% ls() &&
'as.genes.noodles.ranges' %in% ls() &&
class(as.genes.noodles.ranges)=='GRanges' &&	
(
	class(mbd.as.genes.noodles.coverage)=='dgCMatrix' || 
	class(mbd.as.genes.noodles.coverage)=='matrix')
)	noodles.as.genes.tables.loaded <-TRUE
#if they are already in the space, do nothing

if((!noodles.as.genes.tables.loaded) && file.exists('noodles.as.genes.tables.Rda'))
{
	message('Loading as genes noodles...')
	loaded<-load('noodles.as.genes.tables.Rda')
	if (
		'as.genes.noodles.ranges' %in% loaded &&
		'mbd.as.genes.noodles.coverage' %in%  loaded &&
		class(as.genes.noodles.ranges)=='GRanges' &&	
		(
			class(mbd.as.genes.noodles.coverage)=='dgCMatrix' || 
			class(mbd.as.genes.noodles.coverage)=='matrix')
		)	noodles.as.genes.tables.loaded <-TRUE
}
#if we loaded it, we do nothing

if(!noodles.as.genes.tables.loaded)
{
	mbd.bed.dir<-'../../Methylation/bedfiles/'
	chrs<-nucl.chromosomes.hg19()
	mbd.bed.files<-dir(mbd.bed.dir) 
	mbd.bed.files<-mbd.bed.files[grep('All_',mbd.bed.files,invert=TRUE)] # remove two 'All_' files
	mbd.bed.ids<-sapply(strsplit(mbd.bed.files,split='_'),function(x){if(x[2]!='DNA') x[2] else x[3]}) #somhere id in pos 2, somewhere in 3
	mbd.bed.ids<-sub('DNA','',mbd.bed.ids)
	mbd.bed.ids<-sub('PT2','PT',mbd.bed.ids)
	
	noodlen<-10000
	expand<-1 #noodles
	as.genes.noodles.ranges<-unlist(GRangesList(lapply(as.genes.symbols,function(gene) {
		message(gene)
		currentgenename<-gene
		if(gene=='MIR143HG') currentgenename<-'CARMN'
		if(gene=='LRRC48') currentgenename<-'DRC3'
		if(gene=='UTP11L') currentgenename<-'UTP11'
		entrez<-suppressWarnings(select(org.Hs.eg.db,keys=c(currentgenename),columns=c('ENTREZID'),keytype='SYMBOL')$ENTREZID)
		gene.start<-suppressWarnings(min(select(TxDb.Hsapiens.UCSC.hg19.knownGene,keys=c(entrez),keytype='GENEID',columns=c('TXSTART'))$TXSTART))
		gene.end<-suppressWarnings(max(select(TxDb.Hsapiens.UCSC.hg19.knownGene,keys=c(entrez),keytype='GENEID',columns=c('TXEND'))$TXEND))
		gene.chr<-suppressWarnings(select(TxDb.Hsapiens.UCSC.hg19.knownGene,keys=c(entrez),keytype='GENEID',columns=c('TXCHROM'))$TXCHROM)
		noodstart<-gene.start-(gene.start %% noodlen) - expand*noodlen
		noodend<-gene.end-(gene.end %% noodlen)+(expand+1)*noodlen
		gene.noodles<-data.frame('chr'=gene.chr,'start'=seq(noodstart,noodend-noodlen,noodlen),'end'=seq(noodstart+noodlen-1,noodend,noodlen))
		if(currentgenename!=gene) currentgenename<-paste0(currentgenename,'(',gene,')')
		rownames(gene.noodles)<-paste0(currentgenename,'.',gene.noodles$chr,":",gene.noodles$start,'-',gene.noodles$end)
		noodlesGRanges.gene <- GRanges(gene.noodles,seqinfo=nucl.chromosomes.hg19())
	} )))
	
	mbd.as.genes.noodles.coverage<-count.coverage.of.noodles(as.genes.noodles.ranges,paste(mbd.bed.dir,mbd.bed.files,sep='/'),mbd.bed.ids)
	
	save(file='noodles.as.genes.tables.Rda',list=c('as.genes.noodles.ranges','mbd.as.genes.noodles.coverage'))
}

