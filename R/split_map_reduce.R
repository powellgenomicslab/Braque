#' @title Split expression matrix by groups
#' @description Splits an expression matrix from a Seurat object based on
#' a set of columns stored in the metadata
#' @author Jose Alquicira-Hernandez
#' @param object A seurat object
#' @param by A vector providing the names of the columns to be used as groups
#' stored in the \code{@metadata} slot
#' @param slot Expression matrix to split by groups:
#' ("counts", "data", or "scale.data")
#' @return A data.table containing the group variables and the split matrices
#' @export
#' @importFrom Seurat GetAssayData SplitObject
#' @importFrom data.table data.table as.data.table set setnames setattr "%chin%" rbindlist setcolorder
#' @examples
#' data <- Seurat::pbmc_small
#' split_matrix(data, by = c("groups", "RNA_snn_res.1"))

split_matrix <- function(object, by, slot = "counts"){

  if(!is(object, "Seurat"))
    stop("Input `object` must be of 'Seurat' class")


  # Extract gene expression matrix
  m <- GetAssayData(object, slot = slot)

  # Extract metadata and validate "by" variables
  meta.data <- object[[]]
  i <- match(by, names(meta.data))
  if(any(is.na(i))) stop('"', paste0(by[is.na(i)], collapse = ", "),
                         "' not present in metadata")

  barcodes <- rownames(meta.data)
  meta.data <- as.data.table(meta.data, keep.rownames = "barcode")

  # Get number of cells per group
  n_cells <- meta.data[, .N, by]

  # Determine unique values for each column by groups
  meta_vars <- meta.data[, !"barcode"][, .(vars = list(lapply(.SD, unique))), by]

  # Determine if columns can be collapse by group with a single value
  meta_vars <- sapply(meta_vars$vars, sapply, length)
  meta_vars <- apply(meta_vars, 1, function(x) all(x == 1))
  meta_vars <- names(which(meta_vars))

  # Subset meta variables
  meta_vars_cols <- meta_vars
  meta_vars <- meta.data[, c(by, meta_vars), with = FALSE]
  meta_vars <- unique(meta_vars)


  # Subset grouping variables
  meta.data <- meta.data[,c("barcode", by), with = FALSE]

  # Convert grouping variables to character
  for (j in seq_len(ncol(meta.data))){
    if(class(meta.data[[j]]) == 'factor')
      set(meta.data, j = j, value = as.character(meta.data[[j]]))
  }

  set <- list(data = m, meta.data = meta.data)


  # Determine best splitting strategy
  n_groups <- meta.data[, sapply(.SD, function(x) length(unique(x))), .SDcols = by]
  n_groups <- sort(n_groups)

  split_seurat <- function(set, col){
    # Get barcodes per group
    set_meta.data <- split(set$meta.data, set$meta.data[,..col])
    # Get barcodes only
    barcodes <- lapply(set_meta.data, function(x) x[, barcode])

    # Retrive expression values per group
    set_data <- lapply(barcodes,
                       function(x) set$data[, colnames(set$data) %chin% x, drop = FALSE])

    # Create list with expression data and metadata
    res <- mapply(function(x, y){
      list(data = x, meta.data = y)
    }, set_data, set_meta.data, SIMPLIFY = FALSE)

    res

  }

  # Split data by grouping variables

  set <- list(set)

  message("Subseting data...")

  for(var in names(n_groups)){
    message("Splitting by [", var, "]")
    set <- lapply(set, split_seurat, var)
    set <- unlist(set, recursive = FALSE)
  }

  # Gather expression matrices
  set_data <- lapply(set, function(x) x$data)

  # Create metadata/expression matrices data table
  meta_reference <- lapply(set, function(x) unique(x$meta.data[, ..by]))
  meta_reference <- rbindlist(meta_reference)
  meta_reference[, matrix := set_data]

  # Add group meta variables and number of cells
  meta_reference <- merge(meta_reference, meta_vars, by = by)
  meta_reference <- merge(meta_reference, n_cells, by = by)
  setcolorder(meta_reference, c(by, meta_vars_cols, "N", "matrix"))

  # Add grouping information and signature to final object
  #attr(meta_reference, "by") <- meta_order
  setattr(meta_reference, "by", by)

  meta_reference
}


#' @title Map function to each group matrix
#' @description Applies function to each group matrix after running
#' \code{split_matrix}.
#' @author Jose Alquicira-Hernandez
#' @param object A data.table object created with  \code{split_matrix}
#' @param f Function to apply to each expression matrix. Function must be
#' applied to genes (row-wise manner, e.g. Matrix::rowMeans()).
#' @return A data.table containing an extra column called \code{result} containing
#' the output of the provided function
#' @export
#' @importFrom future.apply future_lapply
#' @examples
#' data <- Seurat::pbmc_small
#' group_mats <- split_matrix(data, by = c("groups", "RNA_snn_res.1"))
#' group_mats <- split_matrix(group_mats, Matrix::rowMeans)

map_matrix <- function(object, f){
  object[, result := lapply(matrix, f)]
}


#' @title Gather results across groups
#' @description Combines output results after running \code{map_matrix}. This
#' function assumes the results consist on output vectors generated from
#' applying a function operating in matrices. The grouping hierarchy in the
#' results follows the order provided via \code{by} when running
#' \code{split_matrix}.
#' @author Jose Alquicira-Hernandez
#' @param object A data.table object created and processed with \code{split_matrix}
#' and \code{map_matrix}
#' @return A matrix containing the combined results
#' @export
#' @examples
#' data <- SeuratObject::pbmc_small
#' group_mats <- split_matrix(data, by = c("groups", "RNA_snn_res.1"))
#' group_mats <- map_matrix(group_mats, Matrix::rowMeans)
#' reduce_matrix(group_mats)


reduce_matrix <- function(object, column = TRUE){

  if(!"result" %in% names(object))
    stop("`map_matrix` has not been run yet")

  levels <- attr(object, "by")
  n <- length(levels)

  last_level <- levels[n]

  if(length(levels) == 1){
    levels <- last_level
  }else{
    levels <- levels[-n]
  }

  if(n == 1){
    res <- list(object)
  }else{
    res <- split(object, object[[last_level]])
  }

  res <- lapply(res, function(x){
    m <- as.data.frame(x$result)
    names(m) <- apply(x[, ..levels], 1, paste0, collapse = ".")
    if(!column) m <- t(m)
    m
  })

  res

}


#' @title Gather results across groups
#' @description Combines output results after running \code{map_matrix}. This
#' function assumes the results consist on output vectors generated from
#' applying a function operating in matrices. The grouping hierarchy in the
#' results follows the order provided via \code{by} when running
#' \code{split_matrix}.
#' @author Jose Alquicira-Hernandez
#' @param object A data.table object created and processed with \code{split_matrix}
#' and \code{map_matrix}
#' @return A list containing Seurat objects corresponding to each of the groups
#' @export
#' @examples
#' data <- SeuratObject::pbmc_small
#' group_mats <- split_matrix(data, by = c("groups", "RNA_snn_res.1"))
#' group_mats <- map_matrix(group_mats, Matrix::rowMeans)
#' reduce_seurat(group_mats)


reduce_seurat <- function(object){

  if(!"result" %in% names(object))
    stop("`map_matrix` has not been run yet")

  levels <- attr(object, "by")
  n <- length(levels)

  last_level <- levels[n]

  # Calculate percentage of gene expression per group
  object[, prop_exp := lapply(matrix, function(x) Matrix::rowSums(x > 0) / ncol(x))]

  if(n == 1){
    res <- list(object)
  }else{
    res <- split(object, object[[last_level]])
  }


  res <- lapply(res, function(x){
    m <- as.data.frame(x$result)
    names(m) <- apply(x[, ..levels], 1, paste0, collapse = ".")
    md <- as.data.frame(x[, !c("matrix", "result")])
    rownames(md) <- names(m)
    seurat_obj <- CreateSeuratObject(counts = m, meta.data = md)
    prop_exp <- as.data.frame(x$prop_exp)
    colnames(prop_exp) <- names(m)
    seurat_obj@misc$prop_exp <-  prop_exp
    seurat_obj
  })


}


