fusions<-read.table('../../../WGS/fusions/OnlyInDisease.fullinfo.sorted',stringsAsFactors = FALSE)
the_fusion_samples<-unique(fusions[(fusions[,2]=='MYB' & fusions[,5]=='NFIB') | (fusions[,2]=='NFIB' & fusions[,5]=='MYB'),1])
write(the_fusion_samples,file = 'the.fusion.samples.txt')


