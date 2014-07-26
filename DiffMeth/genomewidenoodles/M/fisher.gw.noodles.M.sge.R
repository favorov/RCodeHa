#the script work without parallel frame if you just run it
#worker-in-sge-array: worker-in-sge-array (gets its worker id from SGE_TASK_ID and the number of workers from SGE_TASK_LAST)
#worker:  worker my-worker-no workers-no (good for any parallel enviroment, including fork)
#combine: combiner workers-no (good for any parallel enviroment, including fork)

args <- commandArgs(trailingOnly = TRUE)

i.am.worker<-FALSE
i.am.combiner<-FALSE

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
	}
	else if ('worker' == args[1] )
	{
		if (length(args)!=3)
			stop('worker is to have two more args')
		if(suppressWarnings(is.na(my.worker.no<-as.integer(args[2]))))
			stop('Second arg is not number')
		if(suppressWarnings(is.na(workers.no<-as.integer(args[3]))))
			stop('Third arg is not number')
		if(my.worker.no<2 || my.worker.no>100)
			stop('My number of worker is a strange number')
		if(workers.no<2 || workers.no>100)
			stop('Number of workers is a strange number')
		if(my.worker.no>workers.no)
			stop('Number of workers is less than my worker no')
		i.am.worker<-TRUE
	}
	else if ('combiner'==args[1])
	{
		if (length(args)!=2)
			stop('combiner is to have one more args')
		if(suppressWarnings(is.na(workers.no<-as.integer(args[2]))))
			stop('Second arg is not number')
		if(workers.no<2 || workers.no>100)
			stop('Number of workers is a strange number')
		i.am.combiner<-TRUE
	}
	else
		stop('First arg is not a \'worker-in-sge-array\',\'worker\' or \'combiner\'')
}

noodles.M.fisher.results.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('noodles.M.fisher.results.Rda'))
	if ('fisher.p.values' %in% load('noodles.M.fisher.results.Rda'))
			noodles.M.fisher.results.loaded<-TRUE
#if we loaded it, we do nothing.

if(!noodles.M.fisher.results.loaded)
{
	if(!sge || i.am.worker)
	{
		#we need this to load
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
		
		
		noodles.M.methylation=noodles.M.methylation[1:60000,] #test
		message('fishering')

		noodles.number<-dim(noodles.M.methylation)[1]

		contrast<-logical(length(bed.ids))
		contrast[grep('HN',bed.ids)]<-TRUE
		if (!sge)
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
		
		fisher.noodles.M.result<-data.frame('fisher.p.values'=numeric(tests.number),'meth.in.normals.ratio'=numeric(tests.number),'meth.in.tumors.ratio'=numeric(tests.number),
			'OR'=numeric(tests.number),'CI_95_L'=numeric(tests.number),'CI_95_H'=numeric(tests.number))

		for (rown in my.worker.start:my.worker.end) 	
		{
			resultrow<-rown+1-my.worker.start
			cotable<-table(as.logical(noodles.M.methylation[rown,]),contrast)
			if(nrow(cotable)==1)#nonmeth
			{
				fisher.noodles.M.result[resultrow,]<-c(1,0,0,NA,NA,NA)
				next
			}
			fisherres<-fisher.test(cotable)
			fisher.noodles.M.result[resultrow,]<-c(fisherres$p.value,cotable[2,2]/cotable[1,2],cotable[2,1]/cotable[1,1],fisherres$estimate,fisherres$conf.int[1],fisherres$conf.int[2])
		}
		message('done\n')

		message('Saving...\n')
		if(!sge)
			save(file='noodles.M.fisher.results.Rda',list=c('fisher.noodles.M.result','tests.number','contrast'))
		else #sge
			save(file=paste('noodles.M.fisher.results.worker.',my.workers.no,'.Rda',sep=''),list=c('fisher.noodles.M.result','noodles.number','tests.number','contrast','my.workers.no','workers.no'))
	}
}

