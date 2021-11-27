#' @title Retrieve proportion of gene expression per sample
#' @description Retrieve proportion of gene expression per sample from a
#' Seurat object created with \code{Braque}
#' @param object \code{Seurat} object
#' @return A data.frame containing proportion of gene expression per sample
#' @importFrom SeuratObject Misc
#' @export
#' @author Jose Alquicira Hernandez

setMethod("prop", signature = "Seurat", definition = function(object){
  Misc(object)$prop
})

