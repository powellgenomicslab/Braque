#' @title Retrieve metadata
#' @description Retrieves sample-level metadata from dge-pair object
#' @param object \code{dge-pair} object
#' @return A data.table containing sample-level metadata
#' @export
#' @author Jose Alquicira Hernandez

setMethod("metadata", signature = "dge-pair", definition = function(object){
  slot(object, name = "metadata")
})

#' @title Retrieve expression data
#' @description Retrieves expression data from dge-pair object
#' @param object \code{dge-pair} object
#' @return A data frame containing expression data
#' @export
#' @author Jose Alquicira Hernandez

setMethod("expr", signature = "dge-pair", definition = function(object){
  slot(object, name = "expr")
})

#' @title Matrix-like accesor for genes and samples
#' @description Matrix-like accesor for genes and samples for dge-pair object
#' @param object \code{dge-pair} object
#' @param i row index (genes)
#' @param j column index (samples)
#' @param drop drop data.frame structure?
#' @return A \code{dge-pair} object with subsetted data
#' @export
#' @author Jose Alquicira Hernandez

setMethod('[', signature(x = "dge-pair"), definition = function(x, i, j, drop = FALSE){
  md <- metadata(x)
  expr <- expr(x)

  expr <- expr[i , j, drop = drop]

  if(!missing(j)){
    if(is.character(j) | is.factor(j)){

      md <- md[sample_condition %chin% j]
      md <- md[match(colnames(expr), md$sample_condition)]
    }else if(is.integer(j) | is.numeric(j) | is.logical(j)){

      md <- md[j, ]
      md <- md[match(colnames(expr), md$sample_condition)]
    }
  }

  x@expr <- expr
  x@metadata <- md

  stopifnot(all(md$sample_condition == colnames(expr)))

  x

})

#' @title Extractor for sample-level columns
#' @description Extractor for sample-level columns for \code{dge-pair} object
#' @param x \code{dge-pair} object
#' @param name Column name
#' @return A vector containing column of interest
#' @export
#' @author Jose Alquicira Hernandez
#'
setMethod('$', signature(x = "dge-pair"), definition = function(x, name){
  md <- metadata(x)
  md[[name]]
})

#' @title Extractor for sample-level columns
#' @description Extractor for sample-level columns for \code{dge-pair} object
#' @param x \code{dge-pair} object
#' @param i Column names or indexes
#' @return A data.table containing column of interest
#' @export
#' @author Jose Alquicira Hernandez

setMethod('[[', signature(x = "dge-pair"), definition = function(x, i){
  md <- metadata(x)
  md[, ..i]
})

#' @title Retrieve gene names
#' @description Retrieve gene names from \code{dge-pair} object
#' @param x \code{dge-pair} object
#' @return A vector containing gene names
#' @export
#' @author Jose Alquicira Hernandez

setMethod('rownames', signature(x = "dge-pair"), definition = function(x){
  expr <- expr(x)
  rownames(expr)
})

#' @title Retrieve sample names
#' @description Retrieve sample names from \code{dge-pair} object
#' @param x \code{dge-pair} object
#' @return A vector containing sample names
#' @export
#' @author Jose Alquicira Hernandez

setMethod('colnames', signature(x = "dge-pair"), definition = function(x){
  expr <- expr(x)
  colnames(expr)
})


#' @title Retrieve sample-level column names
#' @description Retrieve sample-level column namesfrom \code{dge-pair} object
#' @param x \code{dge-pair} object
#' @return A vector containing ssample-level column names
#' @export
#' @author Jose Alquicira Hernandez

setMethod('names', signature(x = "dge-pair"), definition = function(x){
  md <- metadata(x)
  names(md)
})

#' @title Print \code{dge-pair} object
#' @param x \code{dge-pair} object
#' @return Summary information about genes, samples, and contrast/condition
#' @export
#' @author Jose Alquicira Hernandez

setMethod(
  f = 'show',
  signature = 'dge-pair',
  definition = function(object) {
    meta <- metadata(object)
    expr <- expr(object)
    condition <- object@condition_col
    cat(
      "dge-pair object.\n",
      "Groups : ", paste(unique(meta[[condition]]), collapse = ", "), "\n",
      "Genes : ", nrow(expr), "\n",
      "Samples: ", nrow(meta), "\n",
      "Samples per condition: \n\  ", paste(meta[, .N, by = condition][, paste(..condition, N, sep = " = "), .I][[1]], collapse = "\n   "), "\n"
    )
  }
)

#' @title \code{dge-pair} object constructor
#' @param expr Data.frame with rows as genes and columns as samples/individuals
#' representing gene expression data
#' @param metadata A data.table containing sample-level metadata
#' @param sample_col Column name in sample-level metadata containing sample ids
#' @param condition_col Column name in sample-level metadata containing the constrast,
#' for example: cell type population
#' @return \code{dge-pair} object
#' @export
#' @author Jose Alquicira Hernandez

create_dge <- function(expr, meta, sample_col, condition_col){

  if(!inherits(expr, "data.frame")) stop("Expression data must be a data.frame")
  if(!inherits(meta, "data.table")) stop("Metadata must be a data.table object")
  if(!inherits(sample_col, "character")) stop("Sample column name must be a character object")
  if(!inherits(condition_col, "character")) stop("Condition column name must be a character object")

  if(!sample_col %in% names(meta)) stop("Sample column name is not contained in metadata")
  if(!condition_col %in% names(meta)) stop("Condition column name is not contained in metadata")

  if(length(unique(meta[[condition_col]])) != 2) stop("Condition column must contain only two levels")

  sample_condition <- paste(meta[[condition_col]], meta[[sample_col]], sep = ".")
  meta$sample_condition <- sample_condition

  if(!(all(colnames(expr) %in% meta[["sample_condition"]]) & all(meta[["sample_condition"]] %in% colnames(expr))))
    stop("Sample mistmatch between expression data and metadata")

  meta <- meta[match(colnames(expr),  meta[["sample_condition"]]), ]

  stopifnot(all(meta[["sample_condition"]] == colnames(expr)))

  new("dge-pair", metadata = meta, expr = expr, sample_col = sample_col, condition_col = condition_col)

}
