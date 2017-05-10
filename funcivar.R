library(rtracklayer)
library(ggplot2)
library(GenomicRanges)
library(data.table)
library(stringr)


#### FUNCTIONS

#' Read Regions File
#'
#' @param regions.file A path to a file containing SNPs
#' @param search.window Size of window to search vcf for LD calculations, ignored when not doing LD calculations
#'
#' @return
#' @export
#'
#' @examples
ReadRegionsFile <- function(regions.file, search.window = 2e+05) {
  # Reads a tab seperated regions file in the form chr:loc snp_name ethnicity 8:130685457 rs4295627 EUR returns the variables for snp.range, snp.name, and snp.ethno
  snp.regions <- fread(regions.file, data.table = FALSE, showProgress = FALSE)
  snp.region.split <- unlist(strsplit(as.vector(snp.regions[, 1]), ":"))

  snp.chromosome <- as.character(sapply(strsplit(as.vector(snp.regions[, 1]), ":"), function(x) x[1]))
  snp.loc <- as.numeric(sapply(strsplit(as.vector(snp.regions[, 1]), ":"), function(x) x[2]))


  snp.region.start <- round(snp.loc - search.window/2)
  snp.region.end <- round(snp.loc + search.window/2)

  snp.name <- as.character(snp.regions[, 2])
  snp.ethno <- as.character(snp.regions[, 3])
  snp.regions <- data.frame(snp.chromosome, snp.loc, snp.region.start, snp.region.end, snp.name, snp.ethno, stringsAsFactors = FALSE)
  snp.regions$snp.ethno <- toupper(snp.regions$snp.ethno)
  snp.regions$snp.ethno <- str_replace(snp.regions$snp.ethno, "ALL", "AFR,AMR,ASN,EUR")
  snp.regions$snp.ethno <- sapply(str_split(snp.regions$snp.ethno, ","), function(x) {
    x <- paste(unique(x), collapse = ",")
  })
  return(snp.regions)
}



GetBioFeatures <- function(bio.features.loc = NULL, search.term = NULL, biofeatures.genome = "hg19") {
  if (!is.null(bio.features.loc)) {
    bio.features.loc <- str_replace(bio.features.loc, "/$", "")
    peakfiles.pat <- ".bed$|.broadPeak$|.narrowPeak$|.gappedPeak$"
    get.files <- list.files(path = bio.features.loc, pattern = peakfiles.pat, full.names = TRUE, ignore.case = T)
    if (!is.null(search.term)) {
      get.files <- get.files[grepl(search.term, get.files, ignore.case = TRUE)]
    }
    if (length(get.files) >= 1) {
      genome.seqinfo <- Seqinfo(genome = biofeatures.genome)
      bed.list <- lapply(get.files, function(x, s.info) {
        xf <- fread(input = x, header = FALSE, select = c(1:3), skip = "chr", col.names = c("chr", "start", "end"), encoding = "UTF-8", stringsAsFactors = FALSE, data.table = FALSE, showProgress = FALSE)
        name <- sub("(.*)\\..*", "\\1", basename(x))
        xf <- GRanges(seqnames = xf$chr, ranges = IRanges(start = xf$start + 1L, end = xf$end), strand = "*", feature = base::rep.int(name, nrow(xf)), seqinfo = s.info)
        return(xf)
      }, s.info = genome.seqinfo)
      names(bed.list) <- sub("(.*)\\..*", "\\1", basename(get.files))
    } else {
      bed.list <- NULL
      stop("No files matching search term, try a new search term\n", "and make sure that the directory contains files ending in\n", ".bed, .broadPeak, .narrowPeak, or .gappedPeak")
    }
    if (is.null(bed.list)) {
      return(GRangesList())
    } else {
      return(GRangesList(bed.list))
    }
  } else {
    return(NULL)
  }
}

enrich.segments <- function(fg = NULL, bg = NULL,
                            feature = NULL, test = "fisher.test", CI = 0.95,
                            prior = c(a = 1, b = 1), feature.search = NULL,
                            strict.subset = TRUE) {
  ## narrow down features
  if(!is.null(feature.search)) {
    s.feature <- feature[grepl(feature.search, mcols(feature)$state), ]
  } else {
    s.feature <- feature
  }
  ## identify sections not in features
  feature.gaps <- gaps(s.feature)
  newcols <- mcols(s.feature)
  newcols$name <- "no.sample"
  for(column in colnames(newcols)[!(colnames(newcols) %in% "name")]) {
    newcols[, column] <- NA
  }
  trunc.value <- which(c(NROW(newcols), length(feature.gaps)) == min(NROW(newcols), length(feature.gaps)))
  if(trunc.value == 1) {
    trunc.scale <- length(feature.gaps) %% NROW(newcols)
    trunc.value <- length(feature.gaps) / trunc.scale
  } else {
    trunc.value <- length(feature.gaps)
  }
  newcols <- newcols[1:(trunc.value), ]
  mcols(feature.gaps) <- newcols
  ## create enrichment feature
  feature <- c(s.feature, feature.gaps)
  rm(s.feature, feature.gaps)
  samples <- mcols(feature)[, "name"]

  ## Discover Overlaps of FG
  fg.over <- findOverlaps(fg, feature)

  q.fg.over <- queryHits(fg.over)
  s.fg.over <- subjectHits(fg.over)
  fg.over.matrix <- make.overlap.matrix(q.fg.over,
                                        s.fg.over,
                                        samples)
  fg.over.matrix <- fg.over.matrix[, !grepl("no.sample",
                                            colnames(fg.over.matrix))]
  fg.return <- fg
  mcols(fg.return) <- fg.over.matrix

  ## Discover Overlaps of BG
  bg.over <- findOverlaps(bg, feature)
  q.bg.over <- queryHits(bg.over)
  s.bg.over <- subjectHits(bg.over)
  bg.over.matrix <- make.overlap.matrix(q.bg.over,
                                        s.bg.over,
                                        samples)
  bg.over.matrix <- bg.over.matrix[, !grepl("no.sample",
                                            colnames(bg.over.matrix))]

  sample.size <- dim(fg.over.matrix)[[1]]
  sample.stats <- colSums(fg.over.matrix)

  if(is.logical(strict.subset) & strict.subset) {
    total.size <- dim(bg.over.matrix)[[1]]
    total.stats <- colSums(bg.over.matrix)
  } else {
    total.size <- dim(bg.over.matrix)[[1]] + sample.size
    total.stats <- colSums(bg.over.matrix) + sample.stats
  }
  if(test == "fisher.test") {
    enrichment <- data.frame(sample = character(),
                             fg.ratio = numeric(),
                             bg.ratio = numeric(),
                             p.value = numeric(),
                             lower = numeric(),
                             upper = numeric(),
                             q = numeric(),
                             m = numeric(),
                             n = numeric(),
                             k = numeric(),
                             stringsAsFactors = FALSE)

    for(my.sample in seq_along(colnames(bg.over.matrix))) {
      total.success <- total.stats[my.sample]
      sample.success <- sample.stats[my.sample]

      result <- fisher.test(matrix(c(sample.success,
                                     total.success - sample.success,
                                     sample.size - sample.success,
                                     (total.size - total.success) - sample.size + sample.success),
                                   2, 2),
                            conf.level = CI, alternative = "g")

      result <- data.frame(sample = colnames(bg.over.matrix)[my.sample],
                           fg.ratio = sample.success/sample.size,
                           bg.ratio = total.success/total.size,
                           p.value = result$p.value,
                           or = result$estimate,
                           null = result$null.value,
                           lower = result$conf.int[[1]],
                           upper = result$conf.int[[2]],
                           q = sample.success,
                           m = total.success,
                           n = total.size - total.success,
                           k = sample.size,
                           stringsAsFactors = FALSE)
      enrichment <- rbind(enrichment, result)
    }
  } else if(test == "simulation") {
    if(exists("prior") & !is.null(prior)) {
      a <- prior[["a"]]
      b <- prior[["b"]]
    } else {
      stop("no argument present for prior, needed for simulation")
    }
    enrichment <- data.frame(sample = character(),
                             fg.ratio = numeric(),
                             bg.ratio = numeric(),
                             probability = numeric(),
                             difference = numeric(),
                             lower = numeric(),
                             upper = numeric(),
                             fg.success = numeric(),
                             fg.total = numeric(),
                             bg.success = numeric(),
                             bg.total = numeric(),
                             stringsAsFactors = FALSE)
    CI <- (1 - CI)/2
    for(my.sample in seq_along(colnames(bg.over.matrix))) {
      total.success <- total.stats[my.sample]
      sample.success <- sample.stats[my.sample]

      n1 <- sample.size
      y1 <- sample.success
      n2 <- total.size - sample.size
      y2 <- total.success - sample.success
      # SIMULATION
      I = 100000 # simulations
      theta1 = rbeta(I, y1+a, (n1-y1)+b)
      theta2 = rbeta(I, y2+a, (n2-y2)+b)
      diff = theta1-theta2

      # OUTPUT
      quantiles = quantile(diff,c(CI,0.5,1-CI))

      fisher.res <- fisher.test(matrix(c(sample.success,
                                         total.success - sample.success,
                                         sample.size - sample.success,
                                         (total.size - total.success) - sample.size + sample.success),
                                       2, 2),
                                conf.level = CI)

      result <- data.frame(sample = colnames(bg.over.matrix)[my.sample],
                           fg.ratio = sample.success/sample.size,
                           bg.ratio = total.success/total.size,
                           probability = mean(theta1 > theta2),
                           oddsratio = fisher.res$estimate,
                           odds.lower = fisher.res$conf.int[[1]],
                           odds.upper = fisher.res$conf.int[[2]],
                           difference = quantiles[2],
                           lower = quantiles[1],
                           upper = quantiles[3],
                           fg.success = sample.success,
                           fg.total = sample.size,
                           bg.success = total.success,
                           bg.total = total.size,
                           stringsAsFactors = FALSE)
      enrichment <- rbind(enrichment, result)
    }
  } else {
    stop("test not recognized")
  }
  return(list(enrichment = enrichment, overlaps = fg.return))
}

make.overlap.matrix <- function(query.over, subject.over, samples) {
  overlaps <- data.frame(qhits = query.over,
                         samples = samples[subject.over],
                         stringsAsFactors = FALSE)
  non.ol <- unique(samples)[!(unique(samples) %in%
                                samples[subject.over])]
  if(length(non.ol) > 0) {
    overlaps <- rbind(overlaps,
                      data.frame(qhits = NA,
                                 samples = non.ol,
                                 stringsAsFactors = FALSE))
  }
  overlap.table <- table(overlaps)
  overlap.table <- overlap.table[, sort(colnames(overlap.table))]
  output.matrix <- matrix(overlap.table,
                          nrow = nrow(overlap.table),
                          dimnames = list(rownames(overlap.table),
                                          colnames(overlap.table)))
  if(!is(output.matrix, "matrix")) {
    dim(output.matrix) <- c(length(output.matrix), 1)
    dimnames(output.matrix) <- list(rownames(overlap.table),
                                    colnames(overlap.table))
  }
  output.matrix.tmp <- output.matrix > 0
  output.matrix.tmp <- as.integer(output.matrix.tmp)
  output.matrix.tmp <- matrix(output.matrix.tmp, ncol = ncol(output.matrix), nrow = nrow(output.matrix))
  colnames(output.matrix.tmp) <- colnames(output.matrix)
  rownames(output.matrix.tmp) <- rownames(output.matrix)
  return(output.matrix.tmp)
}

plot.beta <- function(result, prior = c(a = 1, b = 1), showprior = FALSE,
                      showdiff = FALSE, ylim = NULL, xlim = NULL,
                      main = "Posterior Distribution of FG and BG") {
  ylim.e <- ylim; xlim.e <- xlim
  x <- seq(0.001, 0.999, length = 100000)
  a1 <- result$fg.success
  b1 <- result$fg.total
  a2 <- result$bg.success
  b2 <- result$bg.total
  fg.beta <- dbeta(x, a1 + prior[["a"]], b1 - a1 + prior[["b"]])
  bg.beta <- dbeta(x, a2 + prior[["a"]], b2 - a2 + prior[["b"]])
  f <- cbind(fg.beta, bg.beta)
  if(showprior) {
    f <- cbind(f, dbeta(x, prior[["a"]], prior[["b"]]))
  }
  if(showdiff) {
    f <- cbind(f, fg.beta - bg.beta)
  }
  if(is.null(ylim.e)) { ylim <- c(0, max(f)) } else { ylim <- ylim.e }
  if(is.null(xlim.e)) { xlim <- c(0, 1) } else { xlim <- xlim.e }
  colors <- c('#66c2a5', '#fc8d62', '#8da0cb')
  matplot(x, f, ylab="density", xlab = "theta", type="l", xlim = xlim, ylim = ylim,
          lty=c(1,1), col = colors[1:ncol(f)], lwd = 2, main = main)
  legend.label.fg <- sprintf("FG Beta(%s, %s)", a1 + prior[["a"]], b1 + prior[["b"]])
  legend.label.bg <- sprintf("BG Beta(%s, %s)", a2 + prior[["a"]], b2 + prior[["b"]])
  legend.label.p <- sprintf("Prior Beta(%s, %s)", prior[["a"]], prior[["b"]])
  legend.label.d <- sprintf("Difference in Beta Distributions")
  legend.label <- c(legend.label.fg, legend.label.bg)
  if(showprior) {
    legend.label <- c(legend.label, legend.label.p)
  }
  if(showdiff) {
    legend.label <- c(legend.label, legend.label.d)
  }
  legend("topright", legend.label,
         col=colors[1:ncol(f)], lty=1, bty = "n")
  invisible(cbind(x, f))
}

change.to.search.genome <- function(granges.object, search.genome) {
  if (Reduce("&", !is.na(genome(granges.object)))) {
    if (identical(genome(granges.object), genome(search.genome))) {
      return(granges.object)
    }
  }
  if (isTRUE(all.equal(seqlevels(granges.object), seqlevels(search.genome)))) {
    seqinfo(granges.object) <- seqinfo(search.genome)
  } else {
    if (seqlevelsStyle(granges.object) != seqlevelsStyle(search.genome)) {
      seqlevelsStyle(granges.object) <- seqlevelsStyle(search.genome)
    }
    normal.xome <- seqlevels(granges.object)[(regexpr("_", seqlevels(granges.object)) < 0)]
    positions <- unlist(sapply(paste0(normal.xome, "$"), grep, seqnames(seqinfo(search.genome))))
    new2oldmap <- rep(NA, length(seqinfo(search.genome)))
    new2oldmap[positions] <- 1:length(positions)
    seqinfo(granges.object, new2old = new2oldmap) <- seqinfo(search.genome)
  }
  return(granges.object)
}
                              
                               
read.SNPsnaps <- function (x, genome="hg19") {
  bkg.snps <- melt(t(fread(x, sep = '\t', header = TRUE, data.table = FALSE)))
  colnames(bkg.snps)<- c("SNP.set", "tag.idx", "match.snp")
  tail (bkg.snps)
  head (bkg.snps)
  odd  <- function(x) x%%2!=0
  even <- function(x) x%%2==0
  bkg.snps$chromosome <- as.character(paste("chr", unlist(strsplit(as.character(bkg.snps$match.snp), split=':'))[odd(1:(length(bkg.snps$SNP.set)*2))], sep=''))
  bkg.snps$snp.chrom  <- as.numeric(unlist(strsplit(as.character(bkg.snps$match.snp), split=':'))[odd(1:(length(bkg.snps$SNP.set)*2))])
  bkg.snps$snp.loc    <- as.numeric(unlist(strsplit(as.character(bkg.snps$match.snp), split=':'))[even(1:(length(bkg.snps$SNP.set)*2))])
  grsnps <- GRanges(seqnames = bkg.snps$chromosome,
                    ranges  = IRanges(start = bkg.snps$snp.loc,
                                      end   = bkg.snps$snp.loc,
                                      names = bkg.snps$match.snp),
                    strand  = '*',
                    setID   = bkg.snps$SNP.set,
                    seqinfo = Seqinfo(genome=genome))
}
