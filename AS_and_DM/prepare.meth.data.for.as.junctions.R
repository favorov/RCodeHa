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
splice.bed.ids.norms<-grep('HN',splice.bed.ids,invert = TRUE)
splice.bed.ids[splice.bed.ids.norms]<-paste0('Normal',splice.bed.ids[splice.bed.ids.norms])
splice.bed.ids<-sub('^[0-9][0-9]','',splice.bed.ids)

list.by.gene<-GRangesList(lapply(junctionRanks,function(junction.gene) {
	colnames(junction.gene)<-splice.bed.ids
	juncRanges.gene <- as.data.frame(strsplit(sub('-',':',rownames(junction.gene)),':'),stringsAsFactors = FALSE)
	juncRanges.gene <- data.frame(t(juncRanges.gene),stringsAsFactors=FALSE)
	rownames(juncRanges.gene)<-rownames(junction.gene)
	colnames(juncRanges.gene)<-c('chr','start','end')
	juncRanges.gene<-cbind(juncRanges.gene,rank=junction.gene)
	juncRanges.gene$start<-as.numeric(juncRanges.gene$start)
	juncRanges.gene$end<-as.numeric(juncRanges.gene$end)
	junctionsGRanges.gene <- GRanges(juncRanges.gene,seqinfo=nucl.chromosomes.hg19())
} ))
#so, we prepared a juction GRanges object... Now, differential.coverage()
mcols(list.by.gene)<-data.frame(gene=names(list.by.gene))
as.junction.list<-unlist(list.by.gene)
mbd.as.junction.coverage<-count.coverage.of.noodles(as.junction.list,paste(mbd.bed.dir,mbd.bed.files,sep='/'),mbd.bed.ids)
