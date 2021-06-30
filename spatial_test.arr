#!/bin/bash
# @ job_name = spat
# @ initialdir = /scratch/devel/mkuhlwilm
# @ output = /home/devel/mkuhlwilm/logs/%j_%a.out
# @ error = /home/devel/mkuhlwilm/logs/%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 2
# @ wall_clock_limit = 10:00:00
# @ memory = 50000
# @ class = normal
# @ requeue = 1
# @ array = 1
##############################################
ID=$SLURM_ARRAY_TASK_ID
module load xz TABIX/0.2.6 gcc/6.3.0 BCFTOOLS/1.6 SAMTOOLS/1.0 BEDTools/2.26.0 bedops;module load gcc/6.3.0;module load GEOS/3.8.1 UDUNITS/2.2.26 R/3.5.0; module load PROJ/7.0.1;module load GDAL/2.4.2;LIBRARY_PATH=/apps/GCC/6.3.0/lib64:$LIBRARY_PATH

infile=/scratch/devel/shan/Bonozee/FASTQs/GATK/all.chr21.g.vcf.gz
labl=test

ip="${infile} ${labl}"
echo $ip
R CMD BATCH --vanilla --slave "--args ${ip}" /home/devel/mkuhlwilm/programs/fulltest.R /home/devel/mkuhlwilm/logs/fulltest_"$labl".log

exit

