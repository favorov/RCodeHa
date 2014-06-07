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

if (!require('doParallel'))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite('doParallel')
	library('doParallel')
}

parallel.workers<-20

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
	source('prepare.gw.noodles.M.R')
}

noodles.M.fisher.results.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('noodles.M.fisher.results.Rda'))
	if ('fisher.p.values' %in% load('noodles.M.fisher.results.Rda'))
			noodles.M.fisher.results.loaded<-TRUE

if(!noodles.M.fisher.results.loaded)
{
	clust<-makeCluster(parallel.workers)
	registerDoParallel(clust)
	message('fishering')

	contrast<-logical(length(bed.ids))
	contrast[grep('HN',bed.ids)]<-TRUE

	tests.number<-length(noodles.M)
	fisher.p.values<-numeric(tests.number)
	meth.in.normals.ratio<-numeric(tests.number)
	meth.in.tumors.ratio<-numeric(tests.number)
	OR<-numeric(tests.number)
	CI_95_L<-numeric(tests.number)
	CI_95_H<-numeric(tests.number)


	foreach (rown = 1:tests.number) %dopar%
	{
		cotable<-table(as.logical(noodles.M.methylation[rown,]),contrast)
		if(nrow(cotable)==1)#nonmeth
		{
			fisher.p.values[rown]<-1.
			meth.in.tumors.ratio[rown]<-0
			meth.in.normals.ratio[rown]<-0
			OR[rown]<-NA
			CI_95_L[rown]<-NA
			CI_95_H[rown]<-NA
		}
		else #calculate
		{
			fisherres<-fisher.test(cotable)
			fisher.p.values[rown]<-fisherres$p.value
			meth.in.tumors.ratio[rown]<-cotable[2,2]/cotable[1,2]
			meth.in.normals.ratio[rown]<-cotable[2,1]/cotable[1,1]
			OR[rown]<-fisherres$estimate
			CI_95_L[rown]<-fisherres$conf.int[1]
			CI_95_H[rown]<-fisherres$conf.int[2]
		}
	}
	stopCluster(clust)
	message('done\n')
	message('Saving...\n')
	save(file='noodles.M.fisher.results.Rda',list=c('fisher.p.values','tests.number','contrast','OR','CI_95_L','CI_95_H'))
}

