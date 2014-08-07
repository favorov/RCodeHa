#/usr/bin/bash
workers=10
for w in 1:10 
do
	echo "Rscript fisher.gw.noodles.M.sge.R worker $workers"
	echo "Rscript fisher.gw.noodles.M.sge.R worker $w $workers" | qsub -cwd -N "fisher.gw.noodles.M.sge.workers_array_$w_$workers"
done
#echo "Rscript qsub fisher.gw.noodles.M.sge.R combiner $workers" > qsub -cwd -hold_jid "fisher.gw.noodles.M.sge.workers_array_$workers" -N "fisher.gw.noodles.M.sge.combiner_$workers"
