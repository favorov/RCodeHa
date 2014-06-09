if (!require('doParallel'))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite('doParallel')
	library('doParallel')
}

parallel.workers<-20

clust<-makeCluster(parallel.workers)
registerDoParallel(clust)
message('testing')

foreach (rown = 1:40) %dopar%
{
	cat('<',file='log.txt',append=TRUE)
	Sys.sleep(1)
	cat('>',file='log.txt',append=TRUE,sep='')
}
stopCluster(clust)
message('done\n')
