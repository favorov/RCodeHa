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

noodles.C.loaded<-FALSE
# we can the whole thing to noodles.1000.Rda
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
	beddir<-'../../../Methylation/bedfiles/'
	noodle.length<-100
	chrs<-nucl.chromosomes.hg19()
	noodles.C<-prepare.uniform.noodles(chrs,noodle.length)
	bedfiles<-dir(beddir) 
	bedfiles<-bedfiles[grep('All_',bedfiles,invert=TRUE)] # remove two 'All_' files
	bed.ids<-sapply(strsplit(bedfiles,split='_'),function(x){if(x[2]!='DNA') x[2] else x[3]}) #somhere id in pos 2, somewhere in 3
	noodles.C.methylation<-CountCoverageOfNoodles(noodles.C,paste0(beddir,bedfiles),bed.ids)
	save(file='noodles.C.Rda',list=c('noodles.C','noodles.C.methylation','bed.ids','noodle.length'))
}

noodles.C.fisher.results.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('noodles.C.fisher.results.Rda'))
	if ('fisher.p.values' %in% load('noodles.C.fisher.results.Rda'))
			noodles.C.fisher.results.loaded<-TRUE

if(!noodles.C.fisher.results.loaded)
{
	message('fishering')

	contrast<-logical(length(bed.ids))
	contrast[grep('HN',bed.ids)]<-TRUE

	tests.number<-length(noodles.C)
	fisher.p.values<-numeric(tests.number)
	meth.in.normals.ratio<-numeric(tests.number)
	meth.in.tumors.ratio<-numeric(tests.number)
	OR<-numeric(tests.number)
	CI_95_L<-numeric(tests.number)
	CI_95_H<-numeric(tests.number)


	for (rown in 1:tests.number)
	{
		cotable<-table(as.logical(noodles.C.methylation[rown,]),contrast)
		if(nrow(cotable)==1)#nonmeth
		{
			fisher.p.values[rown]<-1.
			meth.in.tumors.ratio[rown]<-0
			meth.in.normals.ratio[rown]<-0
			OR[rown]<-NA
			CI_95_L[rown]<-NA
			CI_95_H[rown]<-NA
			next
		}
		fisherres<-fisher.test(cotable)
		fisher.p.values[rown]<-fisherres$p.value
		meth.in.tumors.ratio[rown]<-cotable[2,2]/cotable[1,2]
		meth.in.normals.ratio[rown]<-cotable[2,1]/cotable[1,1]
		OR[rown]<-fisherres$estimate
		CI_95_L[rown]<-fisherres$conf.int[1]
		CI_95_H[rown]<-fisherres$conf.int[2]
	}

	message('done\n')
	message('Saving...\n')
	save(file='noodles.C.fisher.results.Rda',list=c('fisher.p.values','tests.number','contrast','OR','CI_95_L','CI_95_H'))
}

DM.C.noodles.indices<-which(fisher.p.values*tests.number<0.05)

DM.C.noodles<-noodles.C[DM.C.noodles.indices]

DM.C.noodles$p.value<-fisher.p.values[DM.C.noodles.indices]

#DM.C.noodles$OR<-OR[DM.C.noodles.indices]
#DM.C.noodles$CI_95_L<-CI_95_L[DM.C.noodles.indices]
#DM.C.noodles$CI_95_H<-CI_95_H[DM.C.noodles.indices]
DM.C.noodles$ishyper<-CI_95_L[DM.C.noodles.indices]>1
#DM.C.noodles$met.ratio.norm<-meth.in.normals.ratio[DM.C.noodles.indices]
#DM.C.noodles$met.ratio.tumors<-meth.in.tumors.ratio[DM.C.noodles.indices]
DM.C.table<-as(DM.C.noodles,'data.frame')
write.table(DM.C.table,file='DM.C.noodles.bed',quote=FALSE,row.names=FALSE)
