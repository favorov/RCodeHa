#$ -q zappa 
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=6G,h_vmem=12G


. /etc/profile.d/modules.sh

module load sharedapps
module load r/3.3.1
#module load r/3.0.1
