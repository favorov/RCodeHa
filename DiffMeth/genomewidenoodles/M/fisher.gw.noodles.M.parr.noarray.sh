#/usr/bin/bash
workers=50
for ((w=1; w<=$workers; w++)) 
do
	echo "Rscript fisher.gw.noodles.M.sge.R worker $w $workers"
	echo "Rscript fisher.gw.noodles.M.sge.R worker $w $workers" | qsub -cwd -N "fisher.gw.noodles.M.sge.workers.$w.$workers"
done
Rscript qsub fisher.gw.noodles.M.sge.R combiner $workers
