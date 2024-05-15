library(sceasy)
library(reticulate)
library(monocle) 
library(Seurat)
reticulate::use_condaenv('sceasy')

# Argument Handling and Validation
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: convert_seurat_to_anndata.R <input_directory> <output_directory>")
}

input_dir <- args[1]
output_dir <- args[2]

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist:", input_dir)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE) # Create output dir

# Process Each File in the Directory
for (file in list.files(input_dir, pattern = "*.RDS$", full.names = TRUE)) {
  outfile <- file.path(output_dir, sub("\\.RDS$", ".h5ad", basename(file)))
  
  tryCatch({
    cat("Converting:", file, "to", outfile, "\n")
    data <- readRDS(file)  
    sceasy::convertFormat(data, from = "seurat", to = "anndata", outFile = outfile)
  }, error = function(e) {
    cat("Error processing", file, ":", e$message, "\n")
  })
}

cat("Conversion complete.\n")