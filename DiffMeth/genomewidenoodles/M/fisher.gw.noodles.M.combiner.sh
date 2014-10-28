#/usr/bin/bash
workers=50
Rscript fisher.gw.noodles.M.R combiner $workers || exit 1
echo 'cleaning...'
#for ((w=1; w<=$workers; w++)) 
#do
#	rm fisher.gw.noodles.M.sge.workers.$w.$workers.e*?
#	rm fisher.gw.noodles.M.sge.workers.$w.$workers.o*?
#done
