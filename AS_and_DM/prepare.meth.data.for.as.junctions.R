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



mbd.bed.dir<-'../Methylation/bedfiles/'
noodle.length<-1000
chrs<-nucl.chromosomes.hg19()
mbd.bed.files<-dir(mbd.bed.dir) 
mbd.bed.files<-mbd.bed.files[grep('All_',mbd.bed.files,invert=TRUE)] # remove two 'All_' files
mbd.bed.ids<-sapply(strsplit(mbd.bed.files,split='_'),function(x){if(x[2]!='DNA') x[2] else x[3]}) #somhere id in pos 2, somewhere in 3
mbd.bed.ids<-sub('DNA','',mbd.bed.ids)
mbd.bed.ids<-sub('PT2','PT',mbd.bed.ids)

load('RankingJunctions_Sasha.rda')
junction.FLIBM1<-junctionRanks[[1]]
splice.bed.ids<-colnames(junction.FLIBM1)
splice.bed.ids<-sub('_hg19','',splice.bed.ids)
splice.bed.ids<-sub('Normal$','',splice.bed.ids)
splice.bed.ids.norms<-grep('HN',splice.bed.ids,invert = TRUE)
splice.bed.ids[splice.bed.ids.norms]<-paste0('Normal',splice.bed.ids[splice.bed.ids.norms])
splice.bed.ids<-sub('^[0-9][0-9]','',splice.bed.ids)
colnames(junction.FLIBM1)<-splice.bed.ids
juncRanges.FLIBM1 <- as.data.frame(strsplit(sub('-',':',rownames(junction.FLIBM1)),':'),stringsAsFactors = FALSE)
juncRanges.FLIBM1 <- data.frame(t(juncRanges.FLIBM1),stringsAsFactors=FALSE)
rownames(juncRanges.FLIBM1)<-rownames(junction.FLIBM1)
colnames(juncRanges.FLIBM1)<-c('chr','start','end')
#juncRanges.FLIBM1<-cbind(juncRanges.FLIBM1,rank=junction.FLIBM1)
juncRanges.FLIBM1$start<-as.numeric(juncRanges.FLIBM1$start)
juncRanges.FLIBM1$end<-as.numeric(juncRanges.FLIBM1$end)
junctionsGRanges.FLIBM1 <- GRanges(juncRanges.FLIBM1,seqinfo=nucl.chromosomes.hg19())
#so, we prepared a juction GRanges object... Now, differential.coverage()
mbd.as.junction.coverage<-count.coverage.of.noodles(junctionsGRanges.FLIBM1,paste(mbd.bed.dir,mbd.bed.files,sep='/'),mbd.bed.ids)

