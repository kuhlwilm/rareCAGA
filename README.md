# rareCAGA
## rare Chimpanzee Allele Geographic Assignment

This respository contains scripts necessary to understand or replicate the geolocalization strategy described in "Population dynamics and genetic connectivity in recent chimpanzee history", Fontsere et al., Cell Genomics 2022, https://doi.org/10.1016/j.xgen.2022.100133 \
Any use of the scripts, methods and data in this reposity should be referenced by a citation to this publication.

## How to prepare data

If you have good data (e.g. high coverage genome, or capture of chromosome 21), you may process with standard pipelines for genomic data, i.e. mapping with BWA and genotype calling with GATK. Be aware that you need to include the field DP, AD and GQ in the FORMAT column, if these are missing, you may need to annotate them. In principle, rareCAGA works very well on multi-individual VCF files, and merging the data could speed things up.

In case you have low or unknown coverage of a newly sequenced chimpanzee sample, you may perform mapping with BWA (try mapping quality 30 and remove duplicates), and a simple genotye calling with BCFTOOLS. You need a reference genome (human hg19, either the whole assembly or only chr21), then you can obtain working genotypes like this:

```
refgenome=/path/to/refgenome.fa
name=individual_identifier

bcftools mpileup -f ${refgenome} ${name}.bam -r chr21 -a FORMAT/AD,FORMAT/DP -Oz -o ${name}.mpileup.vcf.gz
bcftools call -f GQ -mv -Oz ${name}.mpileup.vcf.gz -o test/${name}.calls.vcf.gz
tabix -f ${name}.calls.vcf.gz

```

This should work with the current code of rareCAGA (September 2023).


## Requirements

You need R (3.5.0 or 4.2.3 worked for me). You need BCFTOOLS (1.6 or higher). For the spatial part, you need GEOS (3.8.1) UDUNITS (2.2.26) PROJ (7.0.1) and GDAL (2.4.2). Also, you need a couple of R libraries to make the spatial inference and plot it to a nice figure. All R packages used are listed below.
```
library("GenomicRanges")
library("gstat")
library("sp")
library("raster")
library("spData")
library("dplyr")
library("sf")
library("tmap")
library("tmaptools")
library("maptools")
library("automap")
library("spatialEco")
library("splancs")
library("ggplot2")
library("rgdal")
library("mapview")
library("TSCS")
```

There are outdated packages, the code will give multiple warning messages, but it should still run (as of September 2023). A future version of rareCAGA will implement different solutions.


## Scripts
Note that only the script *fulltest.R* is relevant for performing the test itself, as well as *spatial_test.arr*, while the other scripts serve as information on how data was processed.\
The scripts are the following:

<b>filtering.R</b>: A VCF file (generated using snpAD) containing genotypes from fecal samples is parsed into R with bcftools, filtered and turned into a simplified genotype matrix of alleles with an acceptable confidence ("medium quality filtering", see the main paper for details).

<b>get_private_sites.R</b>: The genotype matrix is used to determine "locally rare alleles" (see the main paper for detailed descriptions) for 38 sampling locations across the chimpanzee range.

<b>map_preparation.R</b>: Based on the extant chimpanzee range, a grid of possible geographic origin is created, and the geographic coordinates of the 38 sampling locations are also stored in a proper format.

<b>fulltest.R</b>: The actual geolocalization test. A VCF file containing genotypes of one or more test individuals is parsed and filtered, and intersected with the informative rare alleles. Based on the matching rate to each of the 38 sampling locations, kriging is applied for a full spatial distribution across the range, and plotted as output on a map.

<b>spatial_test.arr</b>: Wrapper script to launch the fulltest on a computation cluster.


## Data
The supporting data in the *data* directory are as follows:

<b>colnam6.Robject</b>: R object with metadata on the samples used for the selection of informative sites (only for information on how data was processed)

<b>coln</b>: R object with sample information for the selection of informative sites (only for information on processing)

<b>spatial_stuff</b>: R object with spatial information, including the chimpanzee range from the IUCN database (only for information on processing)

<b>aprivatT9</b>: R object with informative sites, to be used in the *fulltest.R* test script

<b>derived_sites_bed.tar.gz</b>: Archive with bed files containing the positions of informative sites, needed by *fulltest.R*. <b>This needs to be unpacked before using the script!</b>

\
<b>Sequencing</b> data used for this work is available under the study accession number PRJEB46115. Detailed information on processing of raw reads, mapping and genotype calling can be found in the main publication (Fontsere et al., 2022).
