#/usr/bin/bash
workers=50
for ((w=1; w<=$workers; w++)) 
do
	echo "Rscript fisher.gw.noodles.C.R worker $w $workers"
	wrsh=worker.runner.${w}.of.${workers}.sh
	cp script.R.skel.sh ${wrsh}
	echo "Rscript fisher.gw.noodles.C.R worker $w $workers" >> ${wrsh}
	qsub -N "fisher.gw.noodles.C.worker.${w}.of.${workers}" ${wrsh}
done
