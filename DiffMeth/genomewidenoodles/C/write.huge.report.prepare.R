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

if (!suppressWarnings(require('Matrix')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("Matrix")
	library("Matrix")
}

if (!suppressWarnings(require('data.table')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("data.table")
	library("data.table")
}

if(!('all.for.visible.dasha.report.loaded' %in% ls()) || is.na(all.for.visible.dasha.report.loaded)) all.for.visible.dasha.report.loaded<-FALSE
#for quick-develop
if(!all.for.visible.dasha.report.loaded)
{
	message('loading..')
	load('noodles.C.Rda')
	load('noodles.C.fisher.results.Rda')
	all.for.visible.huge.report.loaded<-TRUE
}

rows.no<-dim(fisher.noodles.C.result)[1]
report.interval<-1:rows.no


huge.loaded<-FALSE
if(file.exists('huge.report.frame.Rda'))
	if ('huge.report.frame' %in% load('huge.report.frame.Rda'))
		if (class(huge.report.frame)=='data.frame')
			huge.loaded<-TRUE

if(!huge.loaded)
{
	report.noodles<-noodles.C[report.interval,]
	report.fisher<-fisher.noodles.C.result[report.interval,]
	#actually, it is to develop for little tests.no

	#prepare dataframe
	message('init dataframe')
	report.frame<-data.frame('chr'=as.character(seqnames(report.noodles)),start=start(report.noodles),end=end(report.noodles),stringsAsFactors = FALSE)

	rownames(report.frame)<-paste(report.frame$chr,":",report.frame$start,'-',report.frame$end,sep='')

	message('adding Fisher')
	report.frame<-cbind(report.frame,
			'fisher.p.value'=report.fisher$fisher.p.values,
			'tmr.ratio'=report.fisher$meth.in.tumors.ratio,
			'nor.ratio'=report.fisher$meth.in.normals.ratio,
			'OR'=report.fisher$OR,
			'CI_95_L'=report.fisher$CI_95_L,
			'CI_95_H'=report.fisher$CI_95_H
		)

	message('Looking for closest genes')
	closest.genes<-closest.gene.start.by.interval(report.noodles)

	message('done')


	message('combining')

	report.frame<-cbind(report.frame,elementMetadata(closest.genes)[,c('closest.TSS','pos','dir','dist')])
	message('done\n')

	message('done')

	huge.report.frame<-data.frame(report.frame)

	message('writing')
	
	save(file='huge.report.frame.Rda',list=c('huge.report.frame'))
	
	message('done')
}


