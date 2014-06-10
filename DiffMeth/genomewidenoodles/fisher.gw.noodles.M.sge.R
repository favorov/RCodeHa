#the script work without parrallel if you just run it
#worker: sge worker #workers-no (gets ites worker id from SGE_TASK_ID)
#combine: sge combiner #worker-no

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0)
{
	if(length(args)!=3)
		stop('Argument nomber is not 3')
	if (args[1] != 'sge') 
		stop('First argument is unknown!')
	if(suppressWarnings(is.na(workers.no<-as.integer(args[3]))))
		stop('Third arg is not number')
	if(workers.no<2 || workers.no>100)
		stop('Third arg is a strange number')
	if('worker'==args[2])
	{
		i.am.worker<-TRUE
		if(my.workers.no<-as.integer(Sys.getenv('SGE_TASK_ID', unset = "-1"))<0)
			stop('Worker run not in array')
		if(sge.task.last<-as.integer(Sys.getenv('SGE_TASK_LAST', unset = "-1"))!=workers.no)
			stop('Worker run; SGE_TASK_LAST!=workers.no')
		if(sge.task.last<-as.integer(Sys.getenv('SGE_TASK_FIRST', unset = "-1"))!=1)
			stop('Worker run; SGE_TASK_FIRST!=1')
		if(sge.task.last<-as.integer(Sys.getenv('SGE_TASK_STEP', unset = "-1"))!=1)
			stop('Worker run; SGE_TASK_STEP!=1')
	}
	else
	{
		if('combiner'==args[2])
			i.am.worker<-FALSE
		else
			stop('Second arg is not a \'worker\' or \'combiner\'')
	}
	sge<-TRUE
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

