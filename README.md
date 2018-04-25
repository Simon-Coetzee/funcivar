# funcivar

[![DOI](https://zenodo.org/badge/88302964.svg)](https://zenodo.org/badge/latestdoi/88302964)

example
```{r}
## prepare variants
vcf <- GetVariantsInWindow(file = my.local.vcf, position = my.pos)
```
```{r}
## prepare biofeatures
my.files <- list.files("~/my.biofeatures/", full.names = TRUE)
biofeatures <- GetBioFeatures(files = my.files, genome = "hg19")
## or segmentations
my.segs <- list.files("~/my.segmentations/", full.names = TRUE)
segmentations <- GetSegmentations(files = my.segs, genome = "hg19")
```
```{r}
## set up variants
my.sample.sheet <- read.delim("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel", stringsAsFactors = FALSE)
my.sample.sheet <- my.sample.sheet[, c(1:4)]
test.vcf <- SetPopulation(test.vcf, sample_sheet = my.sample.sheet)
test.vcf <- CalcLD(test.vcf, index = "rs4820988", population = "CEU")
```
```{r}
## find overlaps between variants and features
test.vcf <- GetOverlaps(test.vcf, biofeatures)
ShowOverlaps(test.vcf)
```
```{r}
## split variants into foreground and background based upon LD
test.vcf <- SplitVcfLd(test.vcf,
		       ld = c(metric = "R.squared", cutoff = 0.8, maf = 0.01),
		       strict.subset = TRUE)
```
```{r}
## calculate the enrichment
enrichment <- enrich.variants(test.vcf, biofeatures, feature.type = "biofeatures", return.overlaps = FALSE)
```
```{r}
## plot results
PlotEnrichment(enrichment, value = "difference", block1 = "celltype", color.by = "sample", ncol = 1)
```
