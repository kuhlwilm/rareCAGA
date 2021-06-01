# rareCAGA
## rare Chimpanzee Allele Geographic Assignment

This respository contains scripts necessary to understand or replicate the geolocalization strategy described in Fontsere et al. (to be published and linked).

## Requirements
You need R (3.5.0 worked for me). You need BCFTOOLS (1.6 or higher). For the spatial part, you need GEOS (3.8.1) UDUNITS (2.2.26) PROJ (7.0.1) and GDAL (2.4.2). Also, you need a couple of R libraries to make the spatial inference and plot it to a nice figure. All R packages used are listed below.
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


## Scripts
The scripts are the following:

<b>filtering.R</b>: A VCF file (generated using snpAD) containing genotypes from fecal samples is parsed into R with bcftools, filtered turned into a simplified genotype matrix of "medium-quality" alleles.

<b>get_private_sites.R</b>: The genotype matrix is used to determine "locally rare alleles" (see the main paper for detailed descriptions) for 38 sampling locations across the chimpanzee range.

<b>map_preparation.R</b>: Based on the extant chimpanzee range, a grid of possible geographic origin is created, and the geographic coordinates of the 38 sampling locations are also stored in a proper format.

<b>fulltest.R</b>: The actual geolocalization test. A VCF file containing genotypes of test individuals is parsed and filtered, and intersected with the informative rare alleles. Based on the matching rate to each of the 38 sampling locations, kriging is applied for a full spatial distribution across the range, and plotted on a map.

<b>spatial_test.arr</b>: Wrapper script to launch the fulltest.


## Data
The underlying data can be found here:

(link to FigShare for bed files)

(link to ENA for raw data)



