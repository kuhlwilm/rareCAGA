#!/bin/bash
# @ job_name = prepro
# @ output = ~/logs/%j_%a.out
# @ error = ~/logs/%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 4:00:00
# @ memory = 20000
# @ class = normal
# @ requeue = 1
# @ array = 1

ID=$SLURM_ARRAY_TASK_ID
module load htslib bcftools samtools bwa conda --auto
## need to implement cutadapt as environment, otherwise might be a module, depending on the cluster
conda activate cutadapt-4.4

# define reference genome and sample ID (in array jobs with $ID)
refgenome=/path/to/hg19.fa
name=test

# cut adapter
zcat ${name}.fastq.gz | cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 30 - > ${name}_
cutadapt.fastq

# do mapping to hg19 reference genome
bwa mem -M ${refgenome} "${name}"_cutadapt.fastq | samtools sort | samtools view -b -F 4 -o ${name}.bam
samtools index ${name}.bam
# simple genotype call with bcftools mpileup+call
bcftools mpileup --threads 2 -f ${refgenome} ${name}.bam -r chr21 -a FORMAT/AD,FORMAT/DP -Oz -o ${name}.mpileup.vcf.gz
bcftools call -f GQ -mv -Oz ${name}.mpileup.vcf.gz -o ${name}.calls.vcf.gz
tabix -f ${name}.calls.vcf.gz
# this vcf file can be used as input for rareCAGA

exit

