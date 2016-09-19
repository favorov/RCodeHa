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
	save(file='as.and.dm.tables.Rda',list=c('as.junction.ranges','mbd.as.junction.coverage'))
}
