#!/bin/bash
#
#SBATCH --job-name=prepro
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --time=4:00:00
#SBATCH --output=~/logs/pp_%j.out
#SBATCH --error=~/logs/pp_%j.err
#SBATCH --array=1

ID=$SLURM_ARRAY_TASK_ID
## depending on the system, the different programs might be implemented as modules, otherwise they could be set up as environments
module load htslib bcftools samtools bwa conda --auto
conda activate cutadapt-4.4

# define reference genome and sample ID (in array jobs with $ID and a list)
refgenome=/path/to/hg19.fa
# in case you run several at once, you may use this:
ID=$SLURM_ARRAY_TASK_ID
# file containing the full path to the raw fastq files
## CAUTION: NEED TO ADJUST - this example would use /path/to/name, with corresponding /path/to/name.fastq.gz files
file=$(sed -n ${ID}p /path/to/named_list.txt)
name=$(echo $file | tr "\/" "\n" | tail -1 )

# cut adapter 
if [[ -f ${name}.fastq.gz ]]; then
echo "single end"
zcat ${name}.fastq.gz | cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 30 - > ${name}_cutadapt.fastq.gz
bwa mem -t 2 -M ${refgenome} "${name}"_cutadapt.fastq.gz | samtools sort -@ 2 | samtools view -@ 2 -b -F 4 -o ${name}.bam
fi

if [[ -f ${name}_1.fastq.gz ]]; then
echo "paired end"
cutadapt -a TACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -m 30 -o ${name}_cutadapt_1.fastq.gz -p ${name}_cutadapt_2.fastq.gz ${name}_1.fastq.gz ${name}_2.fastq.gz
bwa mem -t 2 -M ${refgenome} "${name}"_cutadapt_1.fastq.gz "${name}"_cutadapt_2.fastq.gz | samtools sort -@ 2 | samtools view -@ 2 -b -F 4 -o ${name}.bam
fi

# mapping to the hg19 reference genome
samtools index ${name}.bam
# simple genotype calling with bcftools mpileup+call
bcftools mpileup --threads 2 -f ${refgenome} ${name}.bam -r chr21 -a FORMAT/AD,FORMAT/DP -Ou | bcftools call --threads 2 -f GQ -mv -Oz --write-index -o ${name}.calls.vcf.gz
# this vcf file can be used as input for rareCAGA

echo "done"

exit

