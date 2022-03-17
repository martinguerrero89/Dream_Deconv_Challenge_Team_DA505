#!/bin/Rscript
source("cps_functions.R")
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(glmnet))

#### MAIN
option_list <- list(
  make_option(
    "--input",
    action = "store",
    type = "character",
    help = "Set the name of the expresion data file in CSV format"
  ),
  make_option(
    "--finegrain",
    action = "store_true",
    type = "character",
    default = "FALSE",
    help = "Set the type to fine grain"
  ),
  make_option(
    "--coarsegrain",
    action = "store_true",
    type = "character",
    default = "TRUE",
    help = "Set the type to coarse grain [DEFAULT]"
  ),
  make_option(
    "--output",
    action = "store",
    type = "character",
    default = "",
    help = "Set the name of the output data file in CSV format. You can specify a directory name in the ouput and it will be created by the script." 
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  message("[] Parameters missing. Please use --help for look at available parameters.")
  stop()
} else {
  expression_matrix <- read.csv(opt$input, row.names = 1)
  expresion_matrix <- preprocess_em(expresion_matrix)
  ## Create the directory the output will go into
  
  if (opt$coarsegrain == TRUE) {
    predictions <- do_CPS_coarse(expression_matrix)
  } else{
    predictions <- doc_CPS_fine(expression_matrix)
  }
  dir.create(dirname(opt$output), showWarnings = FALSE)
   write.csv(predictions,
              paste0(opt$output),
              row.names = TRUE)
}
