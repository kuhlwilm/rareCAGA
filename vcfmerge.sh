#!/bin/bash
#
#SBATCH --job-name=vcfmerg
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=40GB
#SBATCH --time=8:00:00
#SBATCH --output=~/logs/vm.out
#SBATCH --error=~/logs/vm.err

## if you want to merge a list of VCF files to rareCAGA them all at once, do this
module load bcftools htslib

vlist=/your/vcflist.txt

bcftools merge -L vlist --threads 2 -Oz --write-index -o /output/directory/all.chr21.vcf.gz

exit
