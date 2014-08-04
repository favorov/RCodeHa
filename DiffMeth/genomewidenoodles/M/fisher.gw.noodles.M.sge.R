#the script work without parallel frame if you just run it
#worker-in-sge-array: worker-in-sge-array (gets its worker id from SGE_TASK_ID and the number of workers from SGE_TASK_LAST)
#worker:  worker my-worker-no workers-no (good for any parallel enviroment, including fork)
#combine: combiner workers-no (good for any parallel enviroment, including fork)

noodle.lenght<-1000
noodle.code<-'M'

args <- commandArgs(trailingOnly = TRUE)

i.am.worker<-FALSE
i.am.combiner<-FALSE
i.am.alone<-FALSE

if (length(args) > 0)
{
	if (! args[1] %in% c('worker-in-sge-array','worker','combiner'))
		stop('First argument is unknown!')
	if( 'worker-in-sge-array' == args[1] )
	{
		if (length(args)>1)
			stop('worker-in-sge-array is to be the only arg')
		if(my.worker.no<-as.integer(Sys.getenv('SGE_TASK_ID', unset = "-1"))<0)
			stop('Worker-in-sge-array run not in sge array (no SGE_TASK_ID set)')
		if(task.last<-as.integer(Sys.getenv('SGE_TASK_LAST', unset = "-1"))<0)
			stop('Worker-in-sge-array run not in sge array (no SGE_TASK_LAST set)')
		if(task.first<-as.integer(Sys.getenv('SGE_TASK_FIRST', unset = "-1"))<0)
			stop('Worker-in-sge-array run not in sge array (no SGE_TASK_FIRST set)')
		if(task.first!=1)
			stop('Worker-in-sge-array run; SGE_TASK_FIRST!=1')
		if(task.step<-as.integer(Sys.getenv('SGE_TASK_STEP', unset = "-1"))<0)
			stop('Worker-in-sge-array run not in sge array (no SGE_TASK_STEP set)')
		if(task.step!=1)
			stop('Worker-in-sge-array run; SGE_TASK_STEP!=1')
		i.am.worker<-TRUE
	} else if ('worker' == args[1] )
	{
		if (length(args)!=3)
			stop('worker is to have two more args')
		if(suppressWarnings(is.na(my.worker.no<-as.integer(args[2]))))
			stop('Second arg is not number')
		if(suppressWarnings(is.na(workers.no<-as.integer(args[3]))))
			stop('Third arg is not number')
		if(my.worker.no<1 || my.worker.no>100)
			stop('My number of worker is a strange number')
		if(workers.no<2 || workers.no>100)
			stop('Number of workers is a strange number')
		if(my.worker.no>workers.no)
			stop('Number of workers is less than my worker no')
		i.am.worker<-TRUE
	} else if ('combiner'==args[1])
	{
		if (length(args)!=2)
			stop('combiner is to have one more args')
		if(suppressWarnings(is.na(workers.no<-as.integer(args[2]))))
			stop('Second arg is not number')
		if(workers.no<2 || workers.no>100)
			stop('Number of workers is a strange number')
		i.am.combiner<-TRUE
	} else
		stop('First arg is not a \'worker-in-sge-array\',\'worker\' or \'combiner\'')
} else
	i.am.alone<-TRUE

if (sum (which(c(i.am.alone,i.am.worker,i.am.combiner))) != 1)
	stop ('Something wrong with the self-identification of the script')

if(i.am.alone || i.am.combiner)
{	
	resultfilename<-paste0('noodles.',noodle.code,'.fisher.results.Rda')
	fisher.results.var.name<-paste0('fisher.noodles.',noodle.code,'.result')
} else #i.am.worker
{
	resultfilename<-paste0('noodles.',noodle.code,'.fisher.results.worker.',my.worker.no,'.Rda')
	fisher.results.var.name<-fisher.noodles.result.mat
}

noodles.fisher.results.loaded<-FALSE
# we can load the whole thing


if(file.exists(resultfilename))
	if (fisher.results.var.name %in% load(resultfilename))
			noodles.fisher.results.loaded<-TRUE
#if we loaded it, we do nothing.

if(!noodles.fisher.results.loaded)
{
	if(i.am.alone || i.am.worker) 
	{
		
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

		#we need this to load - we are going to Fisherise
		noodles.loaded<-FALSE
		# we can load the whole thing from noodles.M.Rda

		noodles.file<-paste0('noodles.',noodle.code,'.Rda')
		noodles.methylation.var.name<-paste0('noodles.',noodle.code,'.methylation') # it is a variable name!
		noodles.ranges.var.name<-paste0('noodles.',noodle.code)

		if(file.exists(noodles.file))
		{
			loaded<-load('noodles.M.Rda')
			if (noodles.methylation.var.name %in% loaded) 
				if (class(get(noodles.methylation.var.name))=='data.frame')
					if (noodles.ranges.var.name %in% loaded)
						if(class(get(noodles.ranges.var.name))=='GRanges')
					noodles.loaded<-TRUE
		}
		if(!noodles.loaded)
		{
			if (i.am.worker) 
				stop (paste0('A worker (',my.worker.no,') cannot load noodles+methylation info from ',noodles.methylation.var.name,' and ',noodles.ranges.var.name,'.\n'))
			#if here, i.am.alone
			source(paste0('prepare.gw.noodles.',noodle.code,'.R'))
		}
		
		#assign(noodles.methylation.var.name,get(noodles.methylation.var.name)[1:60000,]) #test
		message('fishering')

		noodles.number<-dim(get(noodles.methylation.var.name))[1]
		#project-specific code; the bed.ids were prepared in 
		contrast<-logical(length(bed.ids))
		contrast[grep('HN',bed.ids)]<-TRUE
		norm.no<-length(which(!contrast))
		tumor.no<-length(which(contrast))

		fishtabs<-as.matrix(prepare.tabulated.fisher(tumor.no,norm.no))

		if (i.am.alone)
		{
			tests.number<-noodles.number
			my.worker.start=1
			my.worker.end=noodles.number
		}
		else
		{
			tests.number<-noodles.number %/% workers.no 
			my.worker.start<-1+tests.number*(my.workers.no-1)
			my.worker.end<-min(my.worker.start+tests.number-1,noodles.number)
		}
	
		message('create result matrix')
		fisher.noodles.result.mat<-matrix(fishtabs[1,],ncol=6,nrow=tests.number,byrow=TRUE)
		
		colnames(fisher.noodles.result.mat)<-c('fisher.p.values','meth.in.normals.ratio','meth.in.tumors.ratio','OR','CI_95_L','CI_95_H') 
		
		revcontrast<-!contrast
		report.every<-tests.number %/% 100
		message('fill result matrix')
		for (rown in 1:tests.number) 	
		{
			therown<-rown+my.worker.start-1
			if (!(rown %% report.every)) message(therown)
			metraw<-get(noodles.methylation.var.name)[therown,]
			aslogic<-as.logical(metraw)
			MY<-sum(aslogic & contrast)
			MN<-sum(aslogic & revcontrast)
			if (0==MN && 0==MY) next
			fishres<-fishtabs[tab.fisher.row.no(tumor.no,norm.no,MY,MN),]
			fisher.noodles.result.mat[therown,]<-fishres
		}
		if (i.am.alone)
		{
			message('converting to dataframe')
			assign(fisher.results.var.name,as.data.frame(fisher.noodles.result.mat))
			message('done\n')
		}
		message('Saving...\n')
		save(file=resultfilename,list=c(fisher.results.var.name,'tests.number','contrast'))
		#in the worker case, 'fisher.noodles.result.mat' == fisher.results.var.name
	}
	else #combiner
	{
		message('Combiner started...\n')
	}
}
