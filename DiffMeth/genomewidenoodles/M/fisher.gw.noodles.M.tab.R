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

	if (!suppressWarnings(require('Differential.Coverage')))
	{
		if (!suppressWarnings(require('devtools')))
		{
			source("http://bioconductor.org/biocLite.R")
			biocLite("devtools")
			library("devtools")
		}
		install_github('Differential.Coverage','favorov')
		#load_all('../../../../../differential.coverage/')
		library('Differential.Coverage')
	}

	
		
	contrast<-logical(length(bed.ids))
	contrast[grep('HN',bed.ids)]<-TRUE

	norm.no<-length(which(!contrast))
	tumor.no<-length(which(contrast))

	fishtabs<-as.matrix(prepare.tabulated.fisher(tumor.no,norm.no))

	#noodles.M.methylation=noodles.M.methylation[1:60000,] #test
	tests.number<-dim(noodles.M.methylation)[1]
	
	#fisher.noodles.M.result<-data.frame('fisher.p.values'=numeric(tests.number),'meth.in.normals.ratio'=numeric(tests.number),'meth.in.tumors.ratio'=numeric(tests.number), 'OR'=numeric(tests.number),'CI_95_L'=numeric(tests.number),'CI_95_H'=numeric(tests.number)) 

	#fisher.noodles.M.result<-data.frame('fisher.p.values'=numeric(0),'meth.in.normals.ratio'=numeric(0),'meth.in.tumors.ratio'=numeric(0), 'OR'=numeric(0),'CI_95_L'=numeric(0),'CI_95_H'=numeric(0)) 
	
	fisher.noodles.M.result.mat<-matrix(ncol=6,nrow=tests.number)
	
	colnames(fisher.noodles.M.result.mat)<-c('fisher.p.values','meth.in.normals.ratio','meth.in.tumors.ratio','OR','CI_95_L','CI_95_H') 
	
	revcontrast<-!contrast
	report.every<-tests.number/1000
	message('create result matrix')
	for (rown in 1:tests.number) 	
	{
		if (!(rown %% report.every)) message(rown)
		theraw<-noodles.M.methylation[rown,]
		aslogic<-as.logical(theraw)
		MY<-sum(aslogic & contrast)
		MN<-sum(aslogic & revcontrast)
		fishres<-fishtabs[norm.no*MY+MN+1,]
		fisher.noodles.M.result.mat[rown,]<-fishres
	}
	message('converting to dataframe')
	fisher.noodles.M.result<-as.data.frame(fisher.noodles.M.result.mat)
	message('done\n')
	message('Saving...\n')
	save(file='noodles.M.fisher.results.Rda',list=c('fisher.noodles.M.result','tests.number','contrast'))
}

