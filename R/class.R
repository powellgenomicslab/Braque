#' @title Definition of 'dge-pair' class
#' @description An S4 class to containing expression data and sample-level metadata
#' @slot expr Data.frame with rows as genes and columns as samples/individuals
#' representing gene expression data
#' @slot metadata A data.table containing sample-level metadata
#' @slot sample_col Column name in sample-level metadata containing sample ids
#' @slot condition_col Column name in sample-level metadata containing the constrast,
#' for example: cell type population
#' @name dge-pair
#' @rdname dge-pair
#' @aliases dge-pair-class
#' @importFrom data.table data.table
#' @exportClass dge-pair

dge_pair <- setClass("dge-pair", representation(
  metadata = "data.table",
  expr = "data.frame",
  sample_col = "character",
  condition_col = "character"),
  prototype(metadata = data.table(),
            expr = data.frame(),
            sample_col = character(),
            condition_col = character())
)
