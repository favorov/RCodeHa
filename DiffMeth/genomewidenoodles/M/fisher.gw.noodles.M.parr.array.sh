#/usr/bin/bash
workers=50
echo "Rscript qsub fisher.gw.noodles.M.sge.R worker-in-sge-array" | qsub -cwd -t 1-$workers -N "fisher.gw.noodles.M.sge.workers.array.$workers"
echo "Rscript qsub fisher.gw.noodles.M.sge.R combiner $workers" | qsub -cwd -hold_jid "fisher.gw.noodles.M.sge.workers_array_$workers" -N "fisher.gw.noodles.M.sge.combiner.$workers"
