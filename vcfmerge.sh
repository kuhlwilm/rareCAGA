#!/bin/bash
#
#SBATCH --job-name=vcfmerg
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=40GB
#SBATCH --time=8:00:00
#SBATCH --output=~/logs/vm.out
#SBATCH --error=~/logs/vm.err

## if you want to merge a list of VCF files to rareCAGA them all at once, do this
module load bcftools htslib

vlist=/your/vcflist.txt

bcftools merge -l vlist --threads 16 -o /output/directory/all.chr21.vcf.gz
tabix -f /output/directory/all.chr21.vcf.gz

exit
