#/usr/bin/bash
workers=50
Rscript fisher.gw.noodles.M.R combiner $workers || exit 1
