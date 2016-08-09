# beadAnalyze
Automated differential expression analysis for Illumina beadarrays (HT12 V4 arrays).

This R script performs automated differential analysis from Illumina's expression  [idat](https://www.bioconductor.org/packages/devel/bioc/vignettes/illuminaio/inst/doc/EncryptedFormat.pdf) files. Default this script assumes idat files are from HT12 V4 arrays, but should with little modification it should work with any arrays.

Make sure you have the following dependencies installed:

[beadarray](https://www.bioconductor.org/packages/release/bioc/html/beadarray.html), [illuminaHumanv4.db](http://bioconductor.org/packages/release/data/annotation/html/illuminaHumanv4.db.html), [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), ggplot2, ggrepel

###Usage:
Clone this repo and source `analyzeBead.R`

Arguments:

`idats` = input IDAT files

`names` = Sample names for each idat file

`condition` = Sample conditions for each idat file

`ref.condition` = reference condition. 

`fdr` = FDR cutoff. Default 0.05

`pltPCA` = If TRUE performs PCA


```{r}
#Example usage with four idat files (2 control and 2 treated).
source("AnalyzeBead.R")
bead.results = beadAnalyze(idats = c("file1.idat","file2.idat","file3.idat","file4.idat"), 
  names = c("control1","control2","treated1","treated2"), 
  condition = c("control","control","treated","treated"), 
  ref.condition = "treated",
  fdr = 0.05, plotPCA = T) 
```
