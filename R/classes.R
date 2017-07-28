#' overlapVCF class
#'
#' @import methods
#' @import VariantAnnotation
#' @export
setClass("overlapVariants",
         contains = "RangedSummarizedExperiment",
         representation = representation(
           genotype = "SimpleList",
           index.snp = "character",
           population = "character",
           overlap.names = "character"))

#' @export
overlapGRanges <- setClass("overlapGRanges",
                           contains = "GRanges",
                           representation = representation(overlap.names = "character"))
#' @export
ldEnrichmentParams <- setClass("ldEnrichmentParams",
                               representation = representation(metric = "character",
                                                               cutoff = "numeric",
                                                               min_maf = "numeric"),

                               validity = function(object) {
                                 errors <- character()
                                 if(!(is.element(object@metric, c("R.squared", "D.prime") & length(object@metric) == 1))) {
                                   msg <- c("metric must be one of R.squared or D.prime")
                                   errors <- c(errors, msg)
                                 }
                                 if(object@cutoff > 1 | object@cutoff < 0) {
                                   msg <- c("cutoff must be a value between 0 and 1")
                                   errors <- c(errors, msg)
                                 }
                                 if(object@min_maf > 1 | object@min_maf < 0) {
                                   msg <- c("cutoff must be a value between 0 and 1")
                                   errors <- c(errors, msg)
                                 }
                                 if (length(errors) == 0) TRUE else errors
                               })

setGeneric("OverlapData", function(object) standardGeneric("OverlapData"))
#' @export
setMethod("OverlapData",
          signature = signature(object = "overlapVCF"),
          function(object) {
            rowRanges(object@VCF)[, object@overlap.names]
          })
#' @export
setMethod("OverlapData",
          signature = signature(object = "overlapGRanges"),
          function(object) {
            rowRanges(object@GRanges)[, object@overlap.names]
          })
# #' @export
# setMethod("show",
#           signature = signature(object = "overlapVCF"),
#           function(object) {
#             cat("overlapVCF object contains:\n")
#             if(length(object@index.snp) > 0){
#               cat("index snp: ", object@index.snp, "\n")
#             } else {
#               cat("index snp: NA\n")
#             }
#             if(length(object@population) > 0)
#               cat("LD has been calculated for variants against index snp: ", object@index.snp)
#             cat("OverlapData(vcf):\n")
#             cat("  GRanges with ",
#                 length(object@overlap.names),
#                 "columns: ",
#                 object@overlap.names, "\n")
#             show(object@VCF)
#           })


