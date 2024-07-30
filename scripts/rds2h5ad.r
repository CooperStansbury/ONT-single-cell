library(sceasy)
library(reticulate)
library(Seurat)
reticulate::use_condaenv('sceasy')

#' Convert a Seurat Object to AnnData Format
#'
#' This script converts a single Seurat (.RDS) object or a list containing a single Seurat object 
#' from an RDS file into an AnnData (.h5ad) object, saving it to the specified output path.
#' It provides enhanced error handling with detailed stack traces, ensuring compatibility
#' with different Seurat versions.
#'
#' @param input_file Path to the Seurat (.RDS) file.
#' @param output_file Path to save the converted AnnData (.h5ad) file.
#'
#' @return Path to the created .h5ad file
#'
convert_seurat_to_anndata <- function(input_file, output_file) {
  # Input Validation
  if (!file.exists(input_file)) {
    stop("Error: Input file does not exist:", input_file)
  }

  # Ensure output directory exists
  output_dir <- dirname(output_file)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Output directory created (or already exists):", output_dir, "\n")

  cat("Converting:", input_file, "to", output_file, "\n")

  tryCatch({
    data <- readRDS(input_file)

    if (inherits(data, "Seurat")) {
      # Single Seurat object
      seurat_obj <- data
    } else if (inherits(data, "list")) {
      # Check for valid Seurat list
      if (length(data) == 1 && inherits(data[[1]], "Seurat")) {
        seurat_obj <- data[[1]]
      } else {
        stop("Error: Input file contains an invalid list. It should contain only one Seurat object.")
      }
    } else {
      stop("Error: Input file does not contain a Seurat object or a valid list of Seurat objects.")
    }

    # Convert Format (No image slot handling needed)
    sceasy::convertFormat(seurat_obj, from = "seurat", to = "anndata", outFile = output_file)

    cat("  Success! Output file created at:", output_file, "\n")
  }, 
  error = function(e) {
    cat("Error:", e$message, "\n")
    cat("Detailed Stack Trace:\n")
    traceback() 
  })
  
  # Return the path to the created h5ad file
  return(output_file)
}

# Command-Line Interface
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  # Improved docstring with examples
  stop(
"Usage: 
  Rscript convert_seurat_to_anndata.R <input_rds_file> <output_h5ad_file>

Example:
  Rscript convert_seurat_to_anndata.R my_data.rds my_data.h5ad"
  )
}

# Execute Conversion
output_path <- convert_seurat_to_anndata(args[1], args[2])
