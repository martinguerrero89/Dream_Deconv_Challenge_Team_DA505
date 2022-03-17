# Dream Deconvolution Challenge Team DA_505

Team DA_505 submissions to DREAM Challenge Tumor Deconvolution

### Installation

1.  Clone the GitHub repo.

<!-- -->

    git clone https://github.com/martinguerrero89/Dream_Deconv_Challenge_Team_DA505.git

1.  Install the following requirements.

<!-- -->

    R >= 4.0.0
    glmnet

### Usage Script Mode:

Move to the folder where you have cloned the repository

    $ cd Dream_Deconv_Challenge_Team_DA505/

You should provide expression data in a comma-separated file with columns associated to sample ID and rows to genes specified using ENSG codes. An example expression matrix is located in the `input` folder.

Then you can run the script `phase1_cps_run.R` for models from **phase1** or `phase2_cps_run.R` for models used during **phase3** and **validation**. The scripts support the following parameters:

    Options:
        --input=INPUT
            Set the name of the expresion data file in CSV format

        --finegrain
            Set the type to fine grain

        --coarsegrain
            Set the type to coarse grain [DEFAULT]

        --output=OUTPUT
            Set the name of the output data file in CSV format. You can specify a directory name in the ouput and it will be created by the script.

        -h, --help
            Show this help message and exit

#### Examples:

For using a model from phase1 with the fine grain you can run the following example in the command line.

    Rscript ./phase1_cps_run.R --input=input/example.csv --output=tmp/coarse_results.csv --finegrain

Alternatively, you can try coarse grain results using the following command:

    Rscript ./phase1_cps_run.R --input=input/example.csv --output=tmp/coarse_results.csv --coarserain

### Using the functions

You can also use the CPS functions in your code. First, you have to load the source for '*Cell Proportion System*' code.

    source("cps_functions.R")

Load expression data from CSV file respecting the format already mentioned.

    expression_matrix <- read.csv("./input/example.csv",row.names = 1)

Run CPS function specifying if 'coarse' or 'fine' grained predictions are desired

    prediction_data <-do_CPS(expression_matrix, type="coarse")

A dataframe with predictions is returned.

|                     | "Sample_1"          | "Sample_2"          |
|---------------------|---------------------|---------------------|
| "B.cells"           | 0.188899122216215   | 0.141777921974995   |
| "CD4.T.cells"       | 0.309831342020004   | 0.318866359701925   |
| "CD8.T.cells"       | 0.0963248441340425  | 0.0712262527901734  |
| "NK.cells"          | -0.0456347415186432 | -0.0517020535308003 |
| "neutrophils"       | 0.160782989477682   | 0.163020596936101   |
| "monocytic.lineage" | 0.0529310969499952  | 0.131204181347804   |
| "fibroblasts"       | 0.188720384375731   | 0.170000233071814   |
| "endothelial.cells" | 0.283468290480592   | 0.270635345222704   |
