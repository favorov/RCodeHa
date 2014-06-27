noodles.M.fisher.results.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('noodles.M.fisher.results.Rda'))
	if ('fisher.p.values' %in% load('noodles.M.fisher.results.Rda'))
			noodles.M.fisher.results.loaded<-TRUE
#if we loaded it, we do nothing

if(!noodles.M.fisher.results.loaded)
{
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

	#if (!suppressWarnings(require('Differential.Coverage')))
	{
		if (!suppressWarnings(require('devtools')))
		{
			source("http://bioconductor.org/biocLite.R")
			biocLite("devtools")
			library("devtools")
		}
		install_github('Differential.Coverage','favorov')
		load_all('../../../../differential.coverage/')
		#library('Differential.Coverage')
	}

	
		
	stop('ok')
	contrast<-logical(length(bed.ids))
	contrast[grep('HN',bed.ids)]<-TRUE

	tests.number<-dim(noodles.M.methylation)[1]

	fisher.noodles.M.result<-data.frame('fisher.p.values'=numeric(tests.number),'meth.in.normals.ratio'=numeric(tests.number),'meth.in.tumors.ratio'=numeric(tests.number),
		'OR'=numeric(tests.number),'CI_95_L'=numeric(tests.number),'CI_95_H'=numeric(tests.number))

	for (rown in 1:tests.number) 	
	{
		cotable<-table(as.logical(noodles.M.methylation[rown,]),contrast)
		if(nrow(cotable)==1)#nonmeth
		{
			fisher.noodles.M.result[rown,]<-c(1,0,0,NA,NA,NA)
			next
		}
		fisherres<-fisher.test(cotable)
		fisher.noodles.M.result[rown,]<-c(fisherres$p.value,cotable[2,2]/cotable[1,2],cotable[2,1]/cotable[1,1],fisherres$estimate,fisherres$conf.int[1],fisherres$conf.int[2])
	}
	message('done\n')
	message('Saving...\n')
	save(file='noodles.M.fisher.results.Rda',list=c('fisher.noodles.M.result','tests.number','contrast'))
}

