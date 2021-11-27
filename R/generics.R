#' @title Retrieve proportion of gene expression per sample
#' @description Retrieve proportion of gene expression per sample from a
#' Seurat object created with \code{Braque}
#' @param object \code{Seurat} object
#' @return A data.frame containing proportion of gene expression per sample
#' @export
#' @author Jose Alquicira Hernandez
#'
setGeneric("prop",  function(object) {
  standardGeneric("prop")
})
