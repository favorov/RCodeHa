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
	#fisher.p.values<-numeric(tests.number)
	#meth.in.normals.ratio<-numeric(tests.number)
	#meth.in.tumors.ratio<-numeric(tests.number)
	#OR<-numeric(tests.number)
	#CI_95_L<-numeric(tests.number)
	#CI_95_H<-numeric(tests.number)
	#here are the names of the fields in the fisher.result dataframe

	#fisheresult<-data.frame('fisher.p.values'=numeric(0),'meth.in.normals.ratio'=numeric(0),'meth.in.tumors.ratio'=numeric(0),
	#			'OR'=numeric(0),'CI_95_L'=numeric(0),'CI_95_H'=numeric(0)

	noodles.M.methylation=noodles.M.methylation[1:60000,]

	fisher.resultee<-foreach (row=iter(noodles.M.methylation, by='row'),.combine=rbind,.multicombine=TRUE) %dopar%
	{
			cotable<-table(as.logical(row),contrast)
			if(nrow(cotable)==1)#nonmeth
				return(c(1,0,0,NA,NA,NA))	
			fisherres<-fisher.test(cotable)
			c(fisherres$p.value,cotable[2,2]/cotable[1,2],cotable[2,1]/cotable[1,1],fisherres$estimate,fisherres$conf.int[1],fisherres$conf.int[2])
	}
	stopCluster(clust)
	message('done\n')
	colnames(fisher.resulte)<-c('fisher.p.values','meth.in.normals.ratio','meth.in.tumors.ratio','OR','CI_95_L','CI_95_H')
	fisher.noodles.M.result<-as(fisher.resulte,'data.frame')
	message('Saving...\n')
	save(file='noodles.M.fisher.results.Rda',list=c('fisher.results','tests.number','contrast'))
}

