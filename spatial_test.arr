#!/bin/bash

#SBATCH --job-name=spatest
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=50GB
#SBATCH --time=2:00:00
#SBATCH --output=~/logs/rarecaga.out
#SBATCH --error=~/logs/rarecaga.err

##############################################
ID=$SLURM_ARRAY_TASK_ID
##module load gcc/6.3.0;module load GEOS/3.8.1 UDUNITS/2.2.26 R/3.5.0; module load PROJ/7.0.1;module load GDAL/2.4.2;LIBRARY_PATH=/apps/GCC/6.3.0/lib64:$LIBRARY_PATH

## as it works on a different cluster
module load R/4.2.3 bcftools conda; conda activate geo ; conda activate gdal;conda activate proj; LD_LIBRARY_PATH=$LD_LIBRARY_PATH

infile=testfile.vcf.gz
labl=test

ip="${infile} ${labl}"
echo $ip
Rscript fulltest.R $ip

echo "done"

exit


