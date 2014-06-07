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

MYF_id <- select(
	org.Hs.eg.db,
	keys=c('MYB'),
	columns=c('ENTREZID'),
	keytype='SYMBOL'
)

genes<- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

stop('qq')

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
	beddir<-'../../../Methylation/bedfiles/'
	noodle.length<-1000
	chrs<-nucl.chromosomes.hg19()
	noodles.M<-prepare.uniform.noodles(chrs,noodle.length)
	bedfiles<-dir(beddir) 
	bedfiles<-bedfiles[grep('All_',bedfiles,invert=TRUE)] # remove two 'All_' files
	bed.ids<-sapply(strsplit(bedfiles,split='_'),function(x){if(x[2]!='DNA') x[2] else x[3]}) #somhere id in pos 2, somewhere in 3
	noodles.M.methylation<-CountCoverageOfNoodles(noodles.M,paste0(beddir,bedfiles),bed.ids)
	save(file='noodles.M.Rda',list=c('noodles.M','noodles.M.methylation','bed.ids','noodle.length'))
}

noodles.M.fisher.results.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('noodles.M.fisher.results.Rda'))
	if ('fisher.p.values' %in% load('noodles.M.fisher.results.Rda'))
			noodles.M.fisher.results.loaded<-TRUE

if(!noodles.M.fisher.results.loaded)
{
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


	for (rown in 1:tests.number)
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
	save(file='noodles.M.fisher.results.Rda',list=c('fisher.p.values','tests.number','contrast','OR','CI_95_L','CI_95_H'))
}

DM.M.noodles.indices<-which(fisher.p.values*tests.number<0.05)

DM.M.noodles<-noodles.M[DM.M.noodles.indices]

DM.M.noodles$p.value<-fisher.p.values[DM.M.noodles.indices]

#DM.M.noodles$OR<-OR[DM.M.noodles.indices]
#DM.M.noodles$CI_95_L<-CI_95_L[DM.M.noodles.indices]
#DM.M.noodles$CI_95_H<-CI_95_H[DM.M.noodles.indices]
DM.M.noodles$ishyper<-CI_95_L[DM.M.noodles.indices]>1
#DM.M.noodles$met.ratio.norm<-meth.in.normals.ratio[DM.M.noodles.indices]
#DM.M.noodles$met.ratio.tumors<-meth.in.tumors.ratio[DM.M.noodles.indices]
DM.M.table<-as(DM.M.noodles,'data.frame')
write.table(DM.M.table,file='DM.M.noodles.bed',quote=FALSE,row.names=FALSE)

con<-file('DM.M.noodles.strict.bed')
export(DM.M.noodles,con,'bed',ignore.strand=TRUE)
close(con)
