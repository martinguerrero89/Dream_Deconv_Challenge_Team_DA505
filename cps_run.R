#!/bin/Rscript
source("cps_functions.R")
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(e1071))

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
    help = "Set the type to fine grain"
  ),
  make_option(
    "--coarsegrain",
    action = "store_true",
    type = "character",
    help = "Set the type to coarse grain [DEFAULT]"
  ),
  make_option(
    "--glmnet",
    action = "store_true",
    type = "character",
    help = "Use glmnet pre-trained model [DEFAULT]"
    ),
  make_option(
    "--svr",
    action = "store_true",
    type = "character",
    help = "Use SVR pre-trained model" ),
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
  message("CPU RUN script.")
  expression_matrix <- read.csv(opt$input, row.names = 1)
  if (!is.null(opt$finegrain)) {
    message("[] Calculating fine grain predictions...")
    if (!is.null(opt$svr)) {
      message(" - using SVR model")
      predictions <- do_CPS_fine_svr(expression_matrix)
    }
    else {
      message(" - using GLMNET model")
      predictions <- do_CPS_fine_glmnet(expression_matrix)
    }
  } else{
    message("[] Calculating coarse grain predictions...")
    if (!is.null(opt$svr )) {
      message(" - using SVR model")
      predictions <- do_CPS_coarse_svr(expression_matrix)
    }
    else{ 
      message(" - using GLMNET model")
      predictions <- do_CPS_coarse_glmnet(expression_matrix)
      }
    }
  message("[] Done.")
  dir.create(dirname(opt$output), showWarnings = FALSE)
   write.csv(predictions,
              paste0(opt$output),
              row.names = TRUE)
  message("[] Predictions saved in ",opt$output)
}
