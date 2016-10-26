if (!suppressWarnings(require('differential.coverage')))
{
	if (!suppressWarnings(require('devtools')))
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("devtools")
		library("devtools")
	}
	install_github('favorov/differential.coverage')
	#load_all('../../../../differential.coverage/')
	library('differential.coverage')
}

noodles.C.loaded<-FALSE
# we can the whole thing to noodles.C.Rda
if(file.exists('noodles.C.Rda'))
{
	loaded<-load('noodles.C.Rda')
	if ('noodles.C.methylation' %in% loaded) 
		if (class(noodles.C.methylation)=='data.frame')
			if ('noodles.C' %in% loaded)
				if(class(noodles.C)=='GRanges')
			noodles.C.loaded<-TRUE
}

if(!noodles.C.loaded)
{
	beddir<-'../../../../Methylation/bedfiles/'
	noodle.length<-100
	chrs<-nucl.chromosomes.hg19()
	noodles.C<-prepare.covering.noodles(chrs,noodle.length)
	bedfiles<-dir(beddir) 
	bedfiles<-bedfiles[grep('All_',bedfiles,invert=TRUE)] # remove two 'All_' files
	bed.ids<-sapply(strsplit(bedfiles,split='_'),function(x){if(x[2]!='DNA') x[2] else x[3]}) #somhere id in pos 2, somewhere in 3
	noodles.C.methylation<-count.coverage.of.noodles(noodles.C,paste0(beddir,bedfiles),bed.ids)
	save(file='noodles.C.Rda',list=c('noodles.C','noodles.C.methylation','bed.ids','noodle.length'))
}
