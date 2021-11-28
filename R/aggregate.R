#' @title Aggregate expression values by groups
#' @description Creates pseudo-bulk gene expression values based on sample groups
#' @author Jose Alquicira-Hernandez
#' @param object A Seurat object
#' @param by A vector providing the names of the columns to be used as groups
#' stored in the \code{@metadata} slot. Last group listed corresponds to the
#' population of interest
#' @param slot Expression matrix to split by groups:
#' ("counts", "data", or "scale.data")
#' @param method Sum or average expression values
#' @return A list containing Seurat objects corresponding to each of the groups
#' @export
#' @importFrom methods is
#' @importFrom SeuratObject GetAssayData CreateSeuratObject CreateAssayObject
#' @importFrom data.table data.table as.data.table set setnames setattr rbindlist setcolorder ":=" transpose
#' @importFrom Matrix sparse.model.matrix Diagonal
#' @examples
#' data <- SeuratObject::pbmc_small
#' aggregate_seurat(data, by = c("groups", "RNA_snn_res.1"))

aggregate_seurat <- function(object, by, slot = "counts", method = c("sum", "mean")){

  # Validate object
  if(!is(object, "Seurat"))
    stop("Input `object` must be of 'Seurat' class")

  # Extract gene expression matrix
  m <- GetAssayData(object, slot = slot)

  # Extract metadata and validate "by" variables
  meta.data <- as.data.table(object[[]], keep.rownames = "barcode")
  i <- match(by, names(meta.data))
  if(any(is.na(i))) stop('"', paste0(by[is.na(i)], collapse = ", "),
                         "' not present in metadata")

  # Validate method
  method <- match.arg(method)

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

  # Create design matrix
  groups <- meta.data[, ..by]
  f <- paste0("~ 0 + ", paste0(names(groups), collapse = ":"))
  mm <- sparse.model.matrix(as.formula(f), data = meta.data, drop.unused.levels = TRUE, sep = "::")

  # Remove group combinations with zero cells
  n_groups <- Matrix::colSums(mm)
  j <- n_groups > 0

  n_groups <- n_groups[j]
  mm <- mm[, j]


  # Get percentage of mean expression
  m_exp <- m > 0
  m_exp <- m_exp %*% mm

  # Get mean weights
  n <- 1 / n_groups
  n <- Diagonal(x = n)
  m_exp <- m_exp %*% n

  # Aggregate data (sum counts)
  m <- m %*% mm

  if(method == "mean"){
    m <- m %*% n
  }

  # Clean up metasample ids
  limit <- length(by) * 3
  i <- seq(from = 3, to = limit, 3)
  md <- as.data.table(transpose(strsplit(colnames(mm), ":")))
  md <- md[, ..i]
  setnames(md, new = by)
  ids <- apply(md, 1, paste0, collapse = ":")
  colnames(m) <- ids -> colnames(m_exp)
  md[, .id := ids]

  # Add index for subsetting
  md[, index := .I]

  # Add number of cells
  n_cells <- Matrix::colSums(mm)
  md[,N := n_cells]

  # Adds sample-level metadata
  md <- merge(md, meta_vars, by = by)

  # Split metasample metadata
  groups <- split(md, md[, by[length(by)], with = FALSE])

  # Subset aggregated data
  m <- lapply(groups, function(x) m[, x$index])
  m_exp <- lapply(groups, function(x) m_exp[, x$index])

  # Format metadata
  groups <- lapply(groups, function(x){
    x[, index := NULL]
    .id <- x$.id
    x <- as.data.frame(x)
    rownames(x) <- .id
    x$.id <- NULL
    x
  })

  # Create Seurat objects
  if(!all(names(m) == names(groups))) stop("ID mismatch")
  m <- mapply(function(exp, meta){
    CreateSeuratObject(exp, meta.data = meta)
  }, m, groups)

  # Add percentage of expression
  if(!all(names(m) == names(m_exp))) stop("ID mismatch")
  m <- mapply(function(x, new_assay){
    x[["prop"]] <- CreateAssayObject(new_assay, key = "prop_")
    x
  }, m, m_exp)

  m

}
