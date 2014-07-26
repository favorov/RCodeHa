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
	
	message('create result matrix')
	fisher.noodles.M.result.mat<-matrix(fishtabs[1,],ncol=6,nrow=tests.number,byrow=TRUE)
	
	colnames(fisher.noodles.M.result.mat)<-c('fisher.p.values','meth.in.normals.ratio','meth.in.tumors.ratio','OR','CI_95_L','CI_95_H') 
	
	revcontrast<-!contrast
	report.every<-tests.number %/% 100
	message('fill result matrix')
	for (rown in 1:tests.number) 	
	{
		if (!(rown %% report.every)) message(rown)
		theraw<-noodles.M.methylation[rown,]
		aslogic<-as.logical(theraw)
		MY<-sum(aslogic & contrast)
		MN<-sum(aslogic & revcontrast)
		if (0==MN && 0==MY) next
		fishres<-fishtabs[tab.fisher.row.no(tumor.no,norm.no,MY,MN),]
		fisher.noodles.M.result.mat[rown,]<-fishres
	}
	message('converting to dataframe')
	fisher.noodles.M.result<-as.data.frame(fisher.noodles.M.result.mat)
	message('done\n')
	message('Saving...\n')
	save(file='noodles.M.fisher.results.Rda',list=c('fisher.noodles.M.result','tests.number','contrast'))
}

