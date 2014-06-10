#/usr/bin/bash
workers=20
workid=qsub -t 1-$workers fisher.dw.noodles.M.sge.run.worker.sh
