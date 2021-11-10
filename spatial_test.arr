#!/bin/bash
# @ job_name = spatest
# @ output = ~/logs/%j_%a.out
# @ error = ~/logs/%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 2
# @ wall_clock_limit = 10:00:00
# @ memory = 50000
# @ class = normal
# @ requeue = 1
# @ array = 1
##############################################
ID=$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0;module load GEOS/3.8.1 UDUNITS/2.2.26 R/3.5.0; module load PROJ/7.0.1;module load GDAL/2.4.2;LIBRARY_PATH=/apps/GCC/6.3.0/lib64:$LIBRARY_PATH

infile=testfile.vcf.gz
labl=test

ip="${infile} ${labl}"
echo $ip
R CMD BATCH --vanilla --slave "--args ${ip}" fulltest.R ~/logs/fulltest_"$labl".log

exit

