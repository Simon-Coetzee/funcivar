#' @description Import variants from VCF or BED files, optionally restricted to a specific window.
#' @param file A character vector describing a file to read
#'
#' @param position An optional GRanges object defining the coordinates of the
#'   file to import
#' @param type A character vector stating that the filetype is either 'vcf' or
#'   'bed'
#' @param genome A character vector describing the genome build against which
#'   the variants were mapped e.g. "hg19", "b37", "hg38", "mm10".
#'
#' @return an object of class VCF or GRanges representing the variants
#' @importFrom VariantAnnotation ScanVcfParam
#' @importFrom Rsamtools TabixFile
#' @importFrom rtracklayer import.bed
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<-
#' @export
GetVariantsInWindow <- function(file, position, genome = "hg19", type = "vcf") {
  if (tolower(type) == "vcf") {
    if(!missing(position)) {
      if(is(position, "GRanges")) {
        params <- ScanVcfParam(which = position)
      } else if(position == "all") {
        params <- ScanVcfParam()
      } else if(is(position, "character")) {
        position <- as(position, "GRanges")
        params <- ScanVcfParam(which = position)
      } else {
        stop("I don't understand the position argument, must be unset, GRanges, or 'all'")
      }
    } else {
      stop("without a position argument, the full vcf will be imported into memory,\n",
           "this can be a very expensive operation. Set the position arg to 'all' to allow")
    }
    vcf <- TabixFile(file)
    variants <- getFILE(file, GetVariantsInWindowVCF, params, genome, N.TRIES = 3L)
    return(variants)
  } else if (tolower(type) == "bed") {
    if(!missing(position) & position != 'all') {
      variants <- import.bed(file, which = position, genome = genome)
    } else {
      variants <- import.bed(file, genome = genome)
    }
    return(variants)
  } else {
    stop("type ", type, " is not yet implemented")
  }
}

#' Set Population in VCF
#' @description Add population specific data to your VCF
#' @param vcf An object of class VCF
#'
#' @param sample_sheet a data.frame with columns describing the sample and
#'   population that you may wish to filter on
#' @return returns the vcf that was passed as an argument with it's colData()
#'   modified to include the population data
#' @importFrom S4Vectors merge
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
#' @export
SetPopulation <- function(vcf, sample_sheet) {
  population <- colData(vcf)
  if (!("Samples" %in% colnames(population))) {
    stop("vcf does not contain sample information")
  }
  population$Samples <- rownames(population)
  sample.col <- which(grepl(pattern = "sample",
                            x = colnames(sample_sheet),
                            ignore.case = TRUE))
  if (length(sample.col) > 0L) {
    sample_sheet <- DataFrame(sample_sheet)
    colnames(sample_sheet)[sample.col] <- "Samples"
    population <- S4Vectors::merge(population, sample_sheet, all.x = TRUE, all.y = FALSE, by = "Samples")
    rownames(population) <- population$Samples
    colData(vcf) <- population
    metadata(vcf)$source <- "funciVar"
    metadata(vcf)$ld <- "none"
    return(vcf)
  } else {
    stop("sample sheet does not contain sample column")
  }
}

#' Calculate LD
#' @description Calculate the Linkage Disequilibrium between an index SNP and the rest of the
#' SNPs in the VCF
#' @param vcf An object of class VCF, optionally processed with SetPopulation()
#'   to allow for calculating LD with a specific population.
#'
#' @param index A character vector that contains a SNP that is present in the
#'   vcf
#' @param population An optional character vector specifying the population upon
#'   which to filter the samples in the vcf
#' @param return A character vector, either 'valid' or 'all'. 'valid' returns
#'   valid SNPs from the vcf only, 'all' attempts to return all variants,
#'   including indels. In the case of 'all', still, only SNPs are used for
#'   calculating LD.
#' @param force Logical indicating that even though LD has already been
#'   calculated upon this VCF, you wish to overwrite that information with the
#'   LD calculation of a new index.
#' @return returns the vcf that was passed as an argument with it's rowRanges()
#'   modified to include the information relative to the index. Including ref,
#'   and alt alleles, allele frequencies (in the population queried), distance
#'   to the index, D prime and R squared.
#' @importFrom VariantAnnotation snpSummary isSNV genotypeToSnpMatrix
#' @importFrom SummarizedExperiment rowRanges rowRanges<-
#' @importFrom snpStats ld
#' @importFrom BiocGenerics start end
#' @import methods
#' @export
CalcLD <- function(vcf, index, population, return = "valid", force = TRUE) {
  ## check input
  if(!is(vcf, "VCF")) {
    stop("vcf object must be of the class VCF")
  }
  if (index %in% rownames(vcf)) {
    if (!isSNV(vcf[rownames(vcf) %in% index, ])) {
      stop("funciVar can only calculate LD for SNVs at this time")
    }
  } else {
    stop("index snp ", index, " not found in current vcf")
  }
  if(missing(population)) {
    population <- "ALL"
    vcf.snv <- isSNV(vcf)
    if (return == "valid") {
      if (any(!vcf.snv)) {
        vcf <- vcf[vcf.snv, ]
        warning("CalcLD only operates on SNVs some of your variants were dropped at this step, set return = 'all', to override this behaviour")
      }
      vcf <- vcf[isSNV(vcf, singleAltOnly = TRUE), ]
    }
  } else {
    samples <- as.data.frame(colData(vcf))
    samples.col <- which(samples == population, arr.ind = TRUE)
    if (nrow(samples.col) > 0L) {
      ## metadata
      samples.col <- samples.col[1, "col"]
      pop.samples <- rownames(samples)[grepl(population, samples[, samples.col])]
      vcf <- vcf[, pop.samples]
      vcf.snv <- isSNV(vcf)
      if (return == "valid") {
        if (any(!vcf.snv)) {
          vcf <- vcf[vcf.snv, ]
          warning("CalcLD only operates on SNVs some of your variants were dropped at this step, set return = 'all', to override this behaviour")
        }
        vcf <- vcf[isSNV(vcf, singleAltOnly = TRUE), ]
      }
    } else {
      stop("your population was not found in your vcf.\nUse setPopulation() to add your sample descriptions, and double check that ",
           population, " is present")
    }
  }
  if(nrow(vcf) >= 1L) {
    vcf.ranges <- rowRanges(vcf)[, c("REF", "ALT")]
    mcols(vcf.ranges) <- c(mcols(vcf.ranges), snpSummary(vcf)[, c("a0Freq", "a1Freq", "HWEpvalue")])
    colnames(mcols(vcf.ranges)) <- c("ref", "alt", "refAlleleFreq", "altAlleleFreq", "HWEpvalue")
    mcols(vcf.ranges)[, "indexSNP"] <- index
    mcols(vcf.ranges)[, "population"] <- population
    mcols(vcf.ranges)[, "distanceToIndex"] <- abs(start(vcf.ranges[index, ]) - start(vcf.ranges))
    ## genotype
    vcf.geno <- genotypeToSnpMatrix(vcf)$genotypes
    vcf.ld <- ld(vcf.geno, vcf.geno[, index], stats = c('D.prime', 'R.squared'))
    mcols(vcf.ranges)$D.prime <- vcf.ld$D.prime[, 1]
    mcols(vcf.ranges)$R.squared <- vcf.ld$R.squared[, 1]
    rowRanges(vcf) <- vcf.ranges
    ## xXx rowData(vcf) <- cbind(rowData(vcf), mcols(vcf.ranges))
    if ("none" %in% metadata(vcf)$ld) {
      metadata(vcf)$ld <- index
    } else if(force) {
      metadata(vcf)$ld <- index
    } else if(!("ld" %in% names(metadata(vcf)))) {
      metadata(vcf)$source <- "funciVar"
      metadata(vcf)$ld <- index
    } else {
      stop("ld has already been calculated on this object for index snp: ", index)
    }
  }
  return(vcf)
}

#' Import Biofeatures
#' @description Imports a set of biofeatures in the form of BED files or derivitive formats
#' for the purpose of annotating variants
#' @param files A character vector containing a list of files to
#'
#' @param genome A character vector describing the genome build against which
#'   the biofeatures were created, e.g. "hg19", "b37", "hg38", "mm10".
#'
#' @return A GRangesList containing the biofeatures imported.
#'
#' @importFrom data.table fread
#' @importFrom stringr str_replace
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<- metadata metadata<- DataFrame
#' @importFrom readr read_tsv cols_only col_character col_integer col_number
#' @importFrom utils head
#' @export
GetBioFeatures <- function(files, genome) {
  if (length(files) < 1L) return(NULL)
  good.files <- sapply(files, function(x) file.exists(x))
  if (sum(good.files) < length(files)) {
    if (sum(!good.files) <= 5L) {
      stop(paste("cannot find the following", sum(!good.files), "files:\n"),
           paste0(files[!good.files], "\n"))
    } else {
      stop(paste("cannot find", sum(!good.files), "files. Some of which are:\n"),
           paste0(head(files[!good.files], n = 5L), "\n"))
    }
  } else {
    genome <- tryCatch(Seqinfo(genome = genome), error = NULL)
    bed.list <- lapply(files,
                       function(file, s.info) {
                         if (grepl(".narrowPeak", file, ignore.case = TRUE)) {
                           col.types <- cols_only(chr = col_character(),
                                                  start = col_integer(),
                                                  end = col_integer(),
                                                  name = col_character(),
                                                  score = col_integer(),
                                                  strand = col_character(),
                                                  signalValue = col_number())
                           col.names <- c("chr", "start", "end", "name", "score", "strand", "signalValue")
                           col.numbers <- c(1:7)
                         } else {
                           col.types <- cols_only(chr = col_character(),
                                                  start = col_integer(),
                                                  end = col_integer())
                           col.names <- c("chr", "start", "end")
                           col.numbers <- c(1:3)
                         }
                         if (any(grepl(".gz$", file))) {
                           xf <- suppressMessages(suppressWarnings(read_tsv(file = file, col_names = col.names,
                                                                            col_types = col.types,
                                                                            progress = FALSE)))
                         } else {
                           xf <- fread(input = file, sep = "\t", header = FALSE,
                                       select = col.numbers, skip = "chr",
                                       col.names = col.names,
                                       encoding = "UTF-8",
                                       stringsAsFactors = FALSE,
                                       data.table = FALSE, showProgress = FALSE)
                         }
                         if (ncol(xf) < 7) xf$signalValue <- NA
                         xf <- try(GRanges(seqnames = xf$chr,
                                           ranges = IRanges(start = xf$start + 1L,
                                                            end = xf$end),
                                           strand = "*",
                                           feature = base::rep.int(basename(file), nrow(xf)),
                                           signalValue = xf$signalValue,
                                           seqinfo = s.info))
                         return(xf)
                       }, s.info = genome)
  }
  if (is.null(bed.list)) {
    return(GRangesList())
  } else {
    bed.list <- GRangesList(bed.list)
    names(bed.list) <- basename(files)
    return(bed.list)
  }
}

#' Imports a set of genome segmentation files, particularly those generated by
#' StatePaintR, for the purpose of annotating variants
#' @param files A character vector containing a list of files to import. These
#'   files need to be BED9 format. see
#'   \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}. And they need to
#'   be completely disjoint intervals.
#'
#' @param genome A character vector describing the genome build against which
#'   the segmentations were created, e.g. "hg19", "b37", "hg38", "mm10".
#'
#' @return A GRangesList containing the segmentations imported.
#' @importFrom rtracklayer import.bed
#' @export
GetSegmentations <- function(files) {
  bed.list <- sapply(files, function(file) {
    bed <- import.bed(file)
  })
  bed.names <- sapply(bed.list, function(bed) {
    bed.name <- tryCatch(bed@trackLine@name, error = NA)
  })
  bed.names[is.na(bed.names)] <- files[is.na(bed.names)]
  bed.list <- Map(format.bed, bed.list, bed.names)
  if (is.null(bed.list)) {
    return(GRangesList())
  } else {
    bed.list <- GRangesList(bed.list)
    names(bed.list) <- bed.names
    return(bed.list)
  }
}

format.bed <- function(bed, name) {
  bed.m <- mcols(bed)
  mcols(bed) <- NULL
  mcols(bed)$sample <- name
  mcols(bed)$state <- bed.m$name
  return(bed)
}

#' Split VCF by Linkage Disequilibrium
#' @description Splits a VCF file into foreground and background variant sets based upon a
#' Linkage Disequilibrium cutoff
#'
#' @param vcf an object of class VCF, particularly one that has been processed
#'   by CalcLD()
#'
#' @param ld a vector containing three fields. These include: \itemize{ \item
#'   metric, either 'R.squared' or 'D.prime' \item cutoff, numeric between 0 and
#'   1, the value of the metric upon which to perform the cutoff \item maf,
#'   numeric between 0 and 1, the minor allele frequency considered to be
#'   included in both foreground and background }
#'
#' @param strict.subset boolean indicating if one wishes the foreground  to be a
#'   strict subset of the background, or an independant set.
#' @return a list of length 2 containing in the 'fg' slot the VCF object of the
#'   foreground SNPs, and the bg slot containing the VCF object of the
#'   background SNPs
#' @export
SplitVcfLd <- function(vcf, ld = c(metric = "R.squared", cutoff = 0.8, maf = 0.01), strict.subset = TRUE) {
  if (!is(vcf, "VCF")) {
    stop("parameter vcf must be a VCF object")
  }
  if (!all(names(ld) %in% c("metric", "cutoff", "maf"))) {
    stop("parameter ld must contain the fields 'metric', 'cutoff', and 'maf'")
  }
  if (!(ld[["metric"]] %in% colnames(mcols(rowRanges(vcf))))) {
    stop("in argument ld$metric: '", ld[['metric']],"' must be present in rowRanges(vcf)")
  }
  if (!(ld[["maf"]] >= 0 && ld[["maf"]] <= 1)) {
    stop("in argument ld$maf: '", ld[['maf']], "' must be between the values of 0 and 1")
  }
  filter <- mcols(rowRanges(vcf))[, "altAlleleFreq"] >= ld[["maf"]]
  fg.filter <- filter & mcols(rowRanges(vcf))[, ld[["metric"]]] >= ld[["cutoff"]]
  if (strict.subset) {
    bg.filter <- fg.filter | mcols(rowRanges(vcf))[, ld[["metric"]]] < ld[["cutoff"]]
  } else {
    bg.filter <- filter & !mcols(rowRanges(vcf))[, ld[["metric"]]] >= ld[["cutoff"]]
  }
  metadata(vcf)$strict.subset <- strict.subset
  bg.filter <- fg.filter | mcols(rowRanges(vcf))[, ld[["metric"]]] < ld[["cutoff"]]
  return(list(fg = vcf[fg.filter & !is.na(fg.filter), ], bg = vcf[bg.filter & !is.na(bg.filter), ]))
}


#' Calculate Enrichment
#' @description Calculates the enrichment of foreground variants vs background
#'   variants against a list of biofeatures.
#' @param variants a list of two GRanges or VCF objects representing the
#'   foregound and background variants to be compared. In the format list(fg =
#'   object, bg = object).
#'
#' @param features A GRangesList of biofeatures or segmentations to calculate
#'   enrichment of the variants against. This field is optional if feature.type
#'   = "biofeatures", and the variants objects have already had SetOverlaps
#'   performed against them.
#' @param feature.type a character vector, either "biofeatures" for generic
#'   feature overlaps, or "segmentations" specifying that the features slot is
#'   populated with disjoint genomic segmentations imported by
#'   GetSegmentations()
#' @param CI numeric, the width of the credible interval for the calculation of
#'   the beta-binomal test for enrichment
#' @param prior numeric vector of length two setting the parameters of the prior
#'   of the beta distribution, in the form, e.g., c(a = 0.5, b = 0.5)
#' @param strict.subset either a boolean value, or "guess" indicating if the
#'   foreground is a strict subset of the background, or an independant set.
#' @param return.overlaps boolean value indicating if, in addition to the
#'   enrichment data.frame, a GRanges object representing the overlaps of the
#'   variants and features be returned.
#'
#' @return either a data.frame indicating the enrichment statistics of the
#'   variant set in the features, or a list including the overlaps of the the
#'   variants and features, in addition to the enrichment data.frame. log1p(oddsratio) are reported from fisher.test
#'
#' @importFrom stats fisher.test formula as.formula quantile rbeta
#' @importFrom GenomicRanges gaps
#' @importMethodsFrom GenomicRanges range
#' @importMethodsFrom IRanges subsetByOverlaps
#' @export
CalculateEnrichment <- function(variants, features, feature.type = "biofeatures",
                            CI = 0.95, prior = c(a = 0.5, b = 0.5),
                            strict.subset = "guess", return.overlaps = FALSE) {
  # set strict.subset
  if (strict.subset == "guess") {
    test.sub <- c(metadata(variants$fg)$strict.subset, metadata(variants$bg)$strict.subset)
    if (all(test.sub, !is.null(test.sub))) {
      strict.subset <- TRUE
    } else if (any(metadata(variants$fg)$strict.subset, metadata(variants$bg)$strict.subset)) {
      stop("foreground and background disagree about whether fg is a strict subset of bg\n",
           "create foreground and background with SplitVcfLD")
    } else if (all(rownames(variants$fg) %in% rownames(variants$bg))) {
      strict.subset <- TRUE
    } else {
      strict.subset <- FALSE
    }
  } else if (!is.logical(strict.subset)) {
    stop("strict.subset must be one of 'guess', TRUE, or FALSE")
  }
  # features
  if (!missing(features)) {
    if (feature.type == "biofeatures") {
      variants$fg <- SetOverlaps(variants$fg, features)
      variants$bg <- SetOverlaps(variants$bg, features)
      fg.features <- colnames(mcols(ShowOverlaps(variants$fg)))
      bg.features <- colnames(mcols(ShowOverlaps(variants$bg)))
      all.features <- union(fg.features, bg.features)
      ## first fg
      if (any(!(is.element(all.features, fg.features)))) {
        if (is(variants$fg, "VCF")) {
          mcols(rowRanges(variants$fg))[, all.features[!(is.element(all.features, fg.features))]] <- 0L
        } else if (is(variants$fg, "GRanges")) {
          mcols(variants$fg)[, all.features[!(is.element(all.features, fg.features))]] <- 0L
        }
      }
      if (is(variants$bg, "VCF")) {
        mcols(rowRanges(variants$fg)) <- cbind(mcols(rowRanges(variants$fg))[, 1:metadata(variants$fg)$overlap.offset-1], mcols(rowRanges(variants$fg))[, all.features])
      } else if (is(variants$fg, "GRanges")) {
        mcols(variants$fg) <- cbind(mcols(variants$fg)[, 1:attributes(variants$fg)$metadata$overlap.offset-1], mcols(variants$fg)[, all.features])
      }
      ## then bg
      if (any(!(is.element(all.features, bg.features)))) {
        if (is(variants$bg, "VCF")) {
          mcols(rowRanges(variants$bg))[, all.features[!(is.element(all.features, bg.features))]] <- 0L
        } else if (is(variants$bg, "GRanges")) {
          mcols(variants$bg)[, all.features[!(is.element(all.features, bg.features))]] <- 0L
        }
      }
      if (is(variants$bg, "VCF")) {
        mcols(rowRanges(variants$bg)) <- cbind(mcols(rowRanges(variants$bg))[, 1:metadata(variants$bg)$overlap.offset-1], mcols(rowRanges(variants$bg))[, all.features])
      } else if (is(variants$bg, "GRanges")) {
        mcols(variants$bg) <- cbind(mcols(variants$bg)[, 1:attributes(variants$bg)$metadata$overlap.offset-1], mcols(variants$bg)[, all.features])
      }
    }
  }
  # feature.type = "biofeatures"
  if (feature.type == "biofeatures") {
    enrichment <- enrich.features(fg = variants$fg, bg = variants$bg, CI = CI, prior = prior, strict.subset = strict.subset)
    if(return.overlaps) {
      overlaps <- GRangesList(foregound.overlaps = rowRanges(variants$fg), background.overlaps = rowRanges(variants$bg))
      return(list(overlaps = overlaps, enrichment = enrichment))
    } else {
      return(enrichment)
    }
  } else if (feature.type == "segmentations") {
    if (missing(features)) {
      stop("include segmentations as the 'features' argument")
    } else {
      if (is(variants$fg, "GRanges") & is(variants$bg, "GRanges")) {
        search.range <- GenomicRanges::union(range(variants$fg), range(variants$bg))
      } else {
        search.range <- GenomicRanges::union(range(rowRanges(variants$fg)), range(rowRanges(variants$bg)))
      }
      if (is(features, "GRangesList")) {
        features <- unlist(features, use.names = FALSE)
      }
      # nfeatures <- gaps(features)
      # nfeatures <- nfeatures[strand(nfeatures) == "*"]
      # mcols(nfeatures)$sample <- "none"
      # mcols(nfeatures)$state <- "unclassified"
      # features <- c(features, nfeatures)
      features <- keepSeqlevels(subsetByOverlaps(features, search.range), seqlevelsInUse(search.range))
      enrichment <- enrich.segments(fg = variants$fg,
                                    bg = variants$bg,
                                    features = features,
                                    CI = CI,
                                    prior = prior,
                                    strict.subset = strict.subset,
                                    return.overlaps = return.overlaps)
    }
  } else {
    stop("feature.type: ", feature.type, "is not availible; try 'biofeatures' or 'segmentations'")
  }
  return(enrichment)
}


#' Assign variants to features
#' @param variants a VCF or GRanges object representing variants to annotate
#'
#' @param features a GRanges or GRangesList object representing features used to
#'   annotate variants.
#'
#' @return returns the variants object that was passed as an argument, VCF
#'   objects have overlaps stored in rowRanges(), GRanges have overlaps stored
#'   in mcols(), ShowOverlaps() can be used to reveal this data.
#' @importFrom S4Vectors to from
#' @importFrom GenomicRanges findOverlaps
#' @importFrom dplyr left_join
#' @export
SetOverlaps <- function(variants, features) {
  nfeatures <- gaps(features)
  nfeatures <- nfeatures[strand(nfeatures) == "*"]
  mcols(nfeatures)$sample <- "none"
  mcols(nfeatures)$state <- "unclassified"
  features <- c(features, nfeatures)
  if (is.null(names(variants))) {
    names(variants) <- as.character(variants)
  }
  overlaps <- findOverlaps(variants, features, ignore.strand = TRUE)
  overlaps <- data.frame(from = names(variants)[from(overlaps)],
                         to = mcols(features)[to(overlaps), "sample"],
                         stringsAsFactors = FALSE)
  resmatrix <- make.overlap.matrix(overlaps)
  resmatrix <- resmatrix[names(variants), !grepl("^none$",
                                                 colnames(resmatrix)), drop = FALSE]
  if(ncol(resmatrix) > 0L) {
    if(is(variants, "GRanges")) {
      overlap.offset <- ncol(mcols(variants)) + 1
      mcols(variants) <- cbind(mcols(variants), DataFrame(resmatrix))
      attributes(variants)$metadata$overlap.offset <- overlap.offset
    } else if (is(variants, "VCF")) {
      overlap.offset <- ncol(mcols(rowRanges(variants))) + 1
      mcols(rowRanges(variants)) <- cbind(mcols(rowRanges(variants)), DataFrame(resmatrix))
      metadata(variants)$overlap.offset <- overlap.offset
    } else {
      stop("variants must be either a GRanges or VCF object, your object is: ", class(variants))
    }
  } else {
    metadata(variants)$overlap.offset <- ncol(mcols(rowRanges(variants))) + 1
    warning("No overlaps were found for this variant/feature combination")
  }
  return(variants)
}

#' Show the overlaps between variants and features
#' @param variants a VCF or GRanges object representing annotated variants
#'
#' @param feature a character list of features about which to display overlaps, this is optional, if set as NULL, all features are returned.
#'
#' @return a GRanges list representing the overlap between variants and features.
#' @export
ShowOverlaps <- function(variants, feature = NULL) {
  if (is(variants, "GRanges")) {
    offset <- attributes(variants)$metadata$overlap.offset
  } else if (is(variants, "VCF")) {
    offset <- metadata(variants)$overlap.offset
    variants <- rowRanges(variants)
  } else {
    stop("variants must be either a VCF or GRanges object")
  }
  if (ncol(mcols(variants)) < offset) {
    return(DataFrame())
  }
  if (is.null(feature)) {
    if (!is.null(offset)) {
      return(variants[, offset:ncol(mcols(variants))])
    } else {
      stop("no overlap data is availible, run SetOverlaps() first")
    }
  } else {
    if (any(feature %in% colnames(mcols(rowRanges(variants))))) {
      return(variants[, feature])
    } else {
      stop("no overlap data is availible for feature search: ", feature)
    }
  }
}

#' Plot Enrichment
#' @description a flexible plotting tool of the enrichment data.frame generated
#'   by CalculateEnrichment()
#' @param variant.enrichment a data.frame of enrichment data generated by
#'   CalculateEnrichment()
#'
#' @param value one of 'difference' or 'oddsratio', thereby plotting the
#'   difference in beta distributions between foreground and background, or the
#'   oddsratio of the comparison.
#' @param block1 character string defining the column of variant.enrichment to
#'   be used as a factor by which to separate the enrichment into blocks
#'   (separate panels) in the plot. Default is NULL, in which case there is no
#'   blocking.
#' @param block2 character string defining the column of variant.enrichment to
#'   be used as a factor by which to separate the enrichment into blocks
#'   (separate panels) in the plot. Default is NULL, in which case there is no
#'   blocking.
#' @param color.by optional character string defining the column of
#'   variant.enrichment to be used to color the values. If this column has hex
#'   color values, in the format #ff0000 they will be used.
#' @param colors optional character vector defining the colors to use for each
#'   value. Must be the same length as the the number of rows in the
#'   variante.enrichment data.frame.
#' @param ncol number of columns to use for facet_wrap if only one block is
#'   defined.
#'
#' @return a ggplot plot object
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr n_distinct
#' @importFrom scales alpha
#' @importFrom graphics plot
#' @importFrom stringr str_length
#' @export
PlotEnrichment <- function(variant.enrichment, value = "difference", block1 = NULL, block2 = NULL, color.by = NULL, colors = NULL, ncol = 1) {

  ## check args
  if (!is(variant.enrichment, "data.frame"))
    stop("variant.enrichment must be a data.frame")
  if (!is.null(block1)) {
    if (!(block1 %in% colnames(variant.enrichment)))
      stop("The block1 argument must either be NULL or a column of variant.enrichment.")
  }
  if (!is.null(block2)) {
    if (!(block2 %in% colnames(variant.enrichment)))
      stop("The block2 argument must either be NULL or a column of variant.enrichment.")
  }
  if (!is.null(color.by)) {
    if (!(color.by %in% colnames(variant.enrichment)))
      stop("The color.by argument must either be NULL or a column of variant.enrichment.")
  }
  if (!(value %in% c("log.odds.ratio", "difference"))) {
    stop("The y argument must be either 'log.odds.ratio' or 'difference'")
  }
  if(value == "log.odds.ratio") value.name <- "log1p(odds ratio)"
  if(value == "difference") value.name <- "difference of foreground \nadbackground distributions"
  if (is.null(color.by)) {
    variant.enrichment$color <- "#ececec"
    variant.enrichment[variant.enrichment$significant, "color"] <- "#ff0000"
    color.is <- "self"
  } else {
    color.vals <- variant.enrichment[, color.by]
    names(color.vals) <- rownames(variant.enrichment)
    color.name <- color.by
    if(!is.null(colors)) {
      if(length(color.vals) == length(colors)) {
        variant.enrichment$color <- colors
        color.is <- "self"
      } else {
        stop("arg 'colors' is not the same length as color.by")
      }
    }
    if (all(grepl("^#", color.vals) & (str_length(color.vals) == 7L))) {
      variant.enrichment$color <- color.vals
      variant.enrichment[!variant.enrichment$significant, "color"] <- alpha(variant.enrichment[!variant.enrichment$significant, "color"], 0.5)
      color.is <- "self"
    } else if (is.character(color.vals) | is.factor(color.vals) | is.logical(color.vals)) {
      color.vals <- as.factor(color.vals)
      if(nlevels(color.vals) < 10L) {
        levels(color.vals) <- brewer.pal(nlevels(color.vals), "Set1")
      } else if(nlevels(color.vals) < 20L) {
        kelly.colours <- c("gray95", "gray13", "gold2", "plum4",
                           "darkorange1", "lightskyblue2", "firebrick",
                           "burlywood3", "gray51", "springgreen4", "lightpink2",
                           "deepskyblue4", "lightsalmon2", "mediumpurple4",
                           "orange", "maroon", "yellow3", "brown4",
                           "yellow4", "sienna4", "chocolate", "gray19")
        levels(color.vals) <- kelly.colours[3:(nlevels(color.vals)+2)]
      } else {
        levels(color.vals) <- colorRampPalette(brewer.pal(9, "Spectral"))(nlevels(color.vals))
      }
      variant.enrichment$color <- as.character(color.vals)
      variant.enrichment[!variant.enrichment$significant, "color"] <- alpha(variant.enrichment[!variant.enrichment$significant, "color"], 0.5)
      color.is <- "self"
    } else {
      variant.enrichment$color <- variant.enrichment[, color.by]
      color.is <- "cont"
    }
  }

  ep <- ggplot(variant.enrichment, aes_string(x = "sample", y = value, group = "color")) + theme_minimal()

  if (n_distinct(variant.enrichment$sample) <= 7L) {
    if(value == "log.odds.ratio") {
      ep <- ep + geom_crossbar(aes(ymin = log.odds.lower, ymax = log.odds.upper, color = color), fatten = 2, width = 0.8)
    } else {
      ep <- ep + geom_crossbar(aes(ymin = lower, ymax = upper, color = color), fatten = 2, width = 0.8)
    }
  } else {
    if(value == "log.odds.ratio") {
      ep <- ep + geom_pointrange(aes(ymin = log.odds.lower, ymax = log.odds.upper, color = color), fatten = 0.5)
    } else {
      ep <- ep + geom_pointrange(aes(ymin = lower, ymax = upper, color = color), fatten = 0.5)
    }
  }
  ep <- ep + scale_x_discrete(name = "Sample")
  if(value == "log.odds.ratio") {
    my.max <- round(max(variant.enrichment$odds.upper), 2) + 0.01
    my.min <- round(min(variant.enrichment$odds.lower), 2) - 0.01
    ep <- ep + scale_y_continuous(name = value.name)
    ep <- ep + geom_hline(yintercept = 0, color = "#c4c4c4", alpha = 0.5) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=log1p(1), alpha=0.1, fill="black")
  } else {
    my.max <- round(max(variant.enrichment$upper), 2) + 0.01
    my.min <- round(min(variant.enrichment$lower), 2) - 0.01
    ep <- ep + scale_y_continuous(name = value.name, breaks = round(seq(from = my.min, to = my.max, length.out = 5), 2))
    ep <- ep + geom_hline(yintercept = 0, color = "#c4c4c4", alpha = 0.5) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, alpha=0.1, fill="black")
  }



  if (all(!is.null(block1), !is.null(block2))) {
    ep <- ep + theme(axis.text.x = element_blank(),
                     legend.position="bottom",
                     panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
                     strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = rel(1.5)),
                     strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = rel(1.5)))
  } else {
    ep <- ep + theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = rel(1.5)),
                     legend.position="bottom",
                     panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
                     strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = rel(1.5)))
  }
  if(!is.null(block1) & !is.null(block2)) {
    my.form <- as.formula(paste(block1, "~", block2))
    ep <- ep + facet_grid(my.form, scales = "free_x", space = "free_x", switch = "x")
  } else if (!is.null(block1) & is.null(block2)) {
    my.form <- as.formula(paste("~", block1))
    ep <- ep + facet_wrap(my.form, ncol = ncol, strip.position = ifelse(ncol > 1L, "top", "right"))
  } else if (is.null(block1) & !is.null(block2)) {
    my.form <- as.formula(paste("~", block2))
    ep <- ep + facet_wrap(my.form, ncol = ncol, strip.position = ifelse(ncol > 1L, "top", "right"))
  }

  if(color.is == "cont") {
    ep <- ep + scale_color_brewer(palette = "Greens")
  } else {
    ep <- ep + scale_color_identity()
  }
  plot(ep)
  return(invisible(ep))
}


###############
### Helpers ###
###############

#' @importFrom GenomeInfoDb seqlevelsStyle keepSeqlevels seqlevelsInUse
#' @importFrom VariantAnnotation readVcf
GetVariantsInWindowVCF <- function(file, param, genome) {
  vcf <- readVcf(file = file, genome = genome, param = param)
  if (nrow(vcf) < 1L) stop("no variants found in interval; \n this is sometimes an error in fetching remote file, try with a local vcf")
  vcf <- keepSeqlevels(vcf, seqlevelsInUse(vcf))
  seqlevelsStyle(vcf) <- "UCSC"
  metadata(vcf)$overlap.offset <- NA
  return(vcf)
}

getFILE <- function(FILE, FUN, ..., N.TRIES=1L) {
  N.TRIES <- as.integer(N.TRIES)
  stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))

  while (N.TRIES > 0L) {
    result <- tryCatch(FUN(FILE, ...), error = identity)
    if (!inherits(result, "error"))
      break
    N.TRIES <- N.TRIES - 1L
  }

  if (N.TRIES == 0L) {
    stop("'getFILE()' failed:",
         "\n  FILE: ", FILE,
         "\n  error: ", conditionMessage(result))
  }

  result
}

make.overlap.matrix <- function(overlaps) {
  overlap.table <- table(overlaps)
  output.matrix <- matrix(overlap.table,
                          nrow = nrow(overlap.table),
                          dimnames = list(rownames(overlap.table),
                                          colnames(overlap.table)))

  if (!is(output.matrix, "matrix")) {
    dim(output.matrix) <- c(length(output.matrix), 1)
    dimnames(output.matrix) <- list(rownames(overlap.table),
                                    colnames(overlap.table))
  }
  output.matrix[output.matrix > 1L] <- 1L
  return(output.matrix)
}

#' @importFrom pbapply pblapply
enrich.segments <- function(fg, bg, features, CI, prior, strict.subset, return.overlaps) {
  if (is(features, "GRangesList")) {
    features <- unlist(features, use.names = FALSE)
  }
  if (all(c("sample", "state") %in% colnames(mcols(features)))) {
    states <- unique(mcols(features)$state)
    if (length(states) > 1L) {
      myapply <- pblapply
    } else {
      myapply <- lapply
    }
    enrichment <- myapply(states, function(state,
                                           local.features = features,
                                           local.fg = fg,
                                           local.bg = bg,
                                           local.CI = CI,
                                           local.prior = prior,
                                           local.strict.subset = strict.subset,
                                           local.return.overlaps = return.overlaps) {
      local.features <- local.features[mcols(local.features)$state %in% state, ]
      local.fg <- SetOverlaps(local.fg, local.features)
      local.bg <- SetOverlaps(local.bg, local.features)
      ## Equalize feature columns
      fg.features <- colnames(mcols(ShowOverlaps(local.fg)))
      bg.features <- colnames(mcols(ShowOverlaps(local.bg)))
      all.features <- union(fg.features, bg.features)
      ## first fg
      if (is(local.fg, "VCF")) {
        if(any(!(is.element(all.features, fg.features)))) {
          mcols(rowRanges(local.fg))[, all.features[!(is.element(all.features, fg.features))]] <- 0L
        }
        mcols(rowRanges(local.fg)) <- cbind(mcols(rowRanges(local.fg))[, 1:metadata(local.fg)$overlap.offset - 1],
                                            mcols(rowRanges(local.fg))[, all.features, drop = FALSE])
      } else if (is(local.fg, "GRanges")) {
        if (any(!(is.element(all.features, fg.features)))) {
          mcols(local.fg)[, all.features[!(is.element(all.features, fg.features))]] <- 0L
        }
        mcols(local.fg) <- cbind(mcols(local.fg)[, 1:attributes(local.fg)$metadata$overlap.offset - 1],
                                 mcols(local.fg)[, all.features, drop = FALSE])
      }
      ## then bg
      if (is(local.bg, "VCF")) {
        if(any(!(is.element(all.features, bg.features)))) {
          mcols(rowRanges(local.bg))[, all.features[!(is.element(all.features, bg.features))]] <- 0L
        }
        mcols(rowRanges(local.bg)) <- cbind(mcols(rowRanges(local.bg))[, 1:metadata(local.bg)$overlap.offset - 1],
                                            mcols(rowRanges(local.bg))[, all.features, drop = FALSE])
      } else if (is(local.bg, "GRanges")) {
        if (any(!(is.element(all.features, bg.features)))) {
          mcols(local.bg)[, all.features[!(is.element(all.features, bg.features))]] <- 0L
        }
        mcols(local.bg) <- cbind(mcols(local.bg)[, 1:attributes(local.bg)$metadata$overlap.offset - 1],
                                 mcols(local.bg)[, all.features, drop = FALSE])
      }
      enrich <- enrich.features(fg = local.fg,
                                bg = local.bg,
                                CI = local.CI,
                                prior = local.prior,
                                strict.subset = local.strict.subset)
      if (local.return.overlaps) {
        return(list(e = enrich, fg.o = ShowOverlaps(local.fg), bg.o = ShowOverlaps(local.bg)))
      } else {
        return(enrich)
      }
    })
    if(return.overlaps) {
      overlaps <- lapply(enrichment, function(x) {
        overlaps <- GRangesList(foregound.overlaps = x$fg.o, background.overlaps = x$bg.o)
        return(overlaps)
      })
      names(overlaps) <- states
      enrichment <- lapply(enrichment, function(x) {
        return(x$e)
      })
    }
    names(enrichment) <- states
    for (state in states) {
      enrichment[[state]]$state <- state
    }
    enrichment <- do.call("rbind", enrichment)
  } else {
    stop("Segmentations must have mcols with fields 'sample' and 'state'")
  }
  if(return.overlaps) {
    return(list(overlaps = overlaps, enrichment = enrichment))
  } else {
    return(enrichment)
  }
}


enrich.features <- function(fg, bg, CI, prior, strict.subset) {
  ## foreground overlaps
  fg.over.matrix <- as.matrix(mcols(ShowOverlaps(fg)))
  rownames(fg.over.matrix) <- rownames(fg)
  ## background overlaps
  bg.over.matrix <- as.matrix(mcols(ShowOverlaps(bg)))
  rownames(bg.over.matrix) <- rownames(bg)

  ## start enrichment
  sample.size <- dim(fg.over.matrix)[[1]]
  sample.stats <- colSums(fg.over.matrix)

  if (strict.subset) {
    total.size <- dim(bg.over.matrix)[[1]]
    total.stats <- colSums(bg.over.matrix)
  } else {
    total.size <- dim(bg.over.matrix)[[1]] + sample.size
    total.stats <- colSums(bg.over.matrix) + sample.stats
  }
  if (exists("prior") & !is.null(prior)) {
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
  for (my.sample in seq_along(colnames(bg.over.matrix))) {
    total.success <- total.stats[my.sample]
    sample.success <- sample.stats[my.sample]

    n1 <- sample.size
    y1 <- sample.success
    n2 <- total.size - sample.size
    y2 <- total.success - sample.success
    # SIMULATION
    I = 100000 # simulations
    theta1 = rbeta(I, y1 + a, (n1 - y1) + b)
    theta2 = rbeta(I, y2 + a, (n2 - y2) + b)
    diff = theta1 - theta2

    # OUTPUT
    quantiles <- try(quantile(diff,c(CI, 0.5, 1 - CI)))
    if(inherits(quantiles, "try-error")) browser()
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
                         log.odds.ratio = log1p(fisher.res$estimate),
                         log.odds.lower = log1p(fisher.res$conf.int[[1]]),
                         log.odds.upper = log1p(fisher.res$conf.int[[2]]),
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
  enrichment$significant <- FALSE
  enrichment[enrichment$probability > (1 - CI) | enrichment$probability < CI, "significant"] <- TRUE
  return(enrichment)
}
