#/usr/bin/bash
workers=50
for ((w=1; w<=$workers; w++)) 
do
	echo "Rscript fisher.gw.noodles.M.R worker $w $workers"
	echo "Rscript fisher.gw.noodles.M.R worker $w $workers" | qsub -cwd -N "fisher.gw.noodles.M.worker.${w}.of.${workers}"
done
