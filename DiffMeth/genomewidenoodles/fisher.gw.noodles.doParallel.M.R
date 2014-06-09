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


	#noodles.M.methylation=noodles.M.methylation[1:60000,]
	#that's why we called it the test
	tests.number<-dim(noodles.M.methylation)[1]

	fisher.resulte<-foreach (row=iter(noodles.M.methylation, by='row')) %dopar%
	{
			cotable<-table(as.logical(row),contrast)
			if(nrow(cotable)==1)#nonmeth
			{
				#return(data.frame('fisher.p.values'=1,'meth.in.normals.ratio'=0,'meth.in.tumors.ratio'=0,'OR'=NA,'CI_95_L'=NA,'CI_95_H'=NA))	
				#fisher.noodles.M.result[rown,]<<-c(1,0,0,NA,NA,NA)
				return(c(1,0,0,NA,NA,NA))
			}
			fisherres<-fisher.test(cotable)
			return(data.frame('fisher.p.values'=fisherres$p.value,'meth.in.normals.ratio'=cotable[2,2]/cotable[1,2],'meth.in.tumors.ratio'=cotable[2,1]/cotable[1,1],'OR'=fisherres$estimate,'CI_95_L'=fisherres$conf.int[1],'CI_95_H'=fisherres$conf.int[2]))	
			#fisher.noodles.M.result[rown,]<<-c(fisherres$p.value,cotable[2,2]/cotable[1,2],cotable[2,1]/cotable[1,1],fisherres$estimate,fisherres$conf.int[1],fisherres$conf.int[2])
			return(c(fisherres$p.value,cotable[2,2]/cotable[1,2],cotable[2,1]/cotable[1,1],fisherres$estimate,fisherres$conf.int[1],fisherres$conf.int[2]))
	}
	stopCluster(clust)
	message('copying')
	fisher.noodles.M.result<-data.frame('fisher.p.values'=numeric(tests.number),'meth.in.normals.ratio'=numeric(tests.number),'meth.in.tumors.ratio'=numeric(tests.number),
		'OR'=numeric(tests.number),'CI_95_L'=numeric(tests.number),'CI_95_H'=numeric(tests.number))

	for(rown in 1:tests.number)
	{
		fisher.noodles.M.result[rown,]=fisher.resulte[[rown]]
	}

	message('Saving...\n')
	save(file='noodles.M.fisher.results.Rda',list=c('fisher.noodles.M.result','tests.number','contrast'))
}

