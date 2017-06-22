# funcivar
example
```{r}
## prepare fg snps
fg.snps <- ReadRegionsFile("fg.probes.txt", search.window = 50)
fg.snps <- GRanges(seqnames = paste0("chr", fg.snps$snp.chromosome),
                   ranges = IRanges(start = fg.snps$snp.loc,
                                    end = fg.snps$snp.loc,
                                    names = fg.snps$snp.name),
                   strand = "*",
                   seqinfo = Seqinfo(genome="hg19"))
fg.snps
```
```{r}
## prepare bg snps
bg.snps <- fread(bgprobes.txt", data.table = F)
colnames(bg.snps) <- c("snp.chromosome", "snp.loc", "blah", "snp.name")
bg.snps <- GRanges(seqnames = paste0("chr", bg.snps$snp.chromosome),
                   ranges = IRanges(start = bg.snps$snp.loc,
                                    end = bg.snps$snp.loc,
                                    names = bg.snps$snp.name),
                   strand = "*",
                   seqinfo = Seqinfo(genome="hg19"))
bg.snps
```
```{r}
## prepare biofeatures
biofeatures <- GetBioFeatures(bio.features.loc = "biofeaturedir/")
## create the correct columns
## mcols must be "name" and "state" - do whatever transform you want to make that work
biofeatures <- lapply(biofeatures, function(x) {name <- mcols(x)$feature; mcols(x) <- NULL; mcols(x)$name <- name; mcols(x)$state <- str_replace(name, ".*_", ""); return(x)})
biofeatures <- do.call(c, unlist(biofeatures, use.names = F))
```
```{r}
## calculate the enrichment
enrich.snps <- enrich.segments(fg = fg.probes,
                               bg = bg.probes,
                               feature = biofeatures,
                               test = "simulation",
                               strict.subset = TRUE)
snp.enrichment <- enrich.snps$enrichment
```
```{r}
## set up significance
snp.enrichment$sig <- "not.significant"
snp.enrichment$sigf <- "not.significant"
snp.enrichment[snp.enrichment$probability > 0.975 | snp.enrichment$probability < 0.025, "sig"] <- "significant.e"
snp.enrichment[snp.enrichment$odds.lower > 1, "sigf"] <- "significant.e"

## set up colors
snp.enrichment$color <- "#ff0000"
snp.enrichment$colorf <- "#ff0000" # fisher tests
snp.enrichment[snp.enrichment$sig == "not.significant", "color"] <- "#ececec"
snp.enrichment[snp.enrichment$sig == "not.significant", "colorf"] <- "#ececec" # fisher tests
```
```{r}
library(ggplot2)
enrich.plot <- ggplot(snp.enrichment, aes(sample, difference, group = color)) +
  geom_crossbar(aes(ymin = lower, ymax = upper, color = color), fatten = 1.2, size = 2) +
  theme_minimal() +
  scale_color_identity() +
  geom_hline(yintercept = 0, color = "#c4c4c4", alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position="none",
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = rel(1.5)),
        strip.text.x = element_text(angle = 270, hjust = 0.5, vjust = 1, size = rel(1.5))) +
  scale_x_discrete(name = "Sample") +
  scale_y_continuous(name = "difference", breaks =  c(0, 0.05, 0.10, 0.15, 0.20))
enrich.plot
```
```{r}
enrich.plot.fisher <- ggplot(snp.enrichment, aes(sample, oddsratio, group = colorf)) +
  geom_crossbar(aes(ymin = odds.lower, ymax = odds.upper, color = colorf), fatten = 1.2, size = 2) +
  theme_minimal() +
  scale_color_identity() +
  geom_hline(yintercept = 0, color = "#c4c4c4", alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position="none",
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = rel(1.5)),
        strip.text.x = element_text(angle = 270, hjust = 0.5, vjust = 1, size = rel(1.5))) +
  geom_hline(yintercept = 1, color = "#ececec", alpha = 0.2) + 
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=1,alpha=0.1,fill="black") +
  scale_x_discrete(name = "Sample") 
enrich.plot.fisher
```