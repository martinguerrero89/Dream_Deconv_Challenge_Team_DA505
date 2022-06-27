# Dream Deconvolution Challenge Team DA_505

Team DA_505 submissions to DREAM Challenge Tumor Deconvolution

### Installation

1.  Clone the GitHub repo.

<!-- -->

    git clone https://github.com/martinguerrero89/Dream_Deconv_Challenge_Team_DA505.git

<!-- -->

2.  Install the following requirements.

    R >= 4.0.0
    glmnet
    e1071
    optparse #if you want to use the scripts

3. Installation time

The expected installation time on a AMD Ryzen 7 with 32Gb running Linux OS is

 1m (install base R) + 3m (packages) + github

### Usage Script Mode:

Move to the folder where you have cloned the repository

    $ cd Dream_Deconv_Challenge_Team_DA505/

You should provide expression data in a comma-separated file with columns associated to sample ID and rows to genes specified using ENSG codes. An example expression matrix is located in the `input` folder.

Then you can run the script `cps_run.R` The scripts support the following parameters:

    Options:
        --input=INPUT
            Set the name of the expresion data file in CSV format

        --finegrain
            Set the type to fine grain

        --coarsegrain
            Set the type to coarse grain

        --glmnet
            Use glmnet pre-trained model

        --svr
            Use SVR pre-trained model

        --output=OUTPUT
            Set the name of the output data file in CSV format. You can specify a directory name in the ouput and it will be created by the script.

        -h, --help
            Show this help message and exit

#### Examples:

For using the GLMNET model (default) for fine grain results you can run the following example in the command line.

    Rscript ./_cps_run.R --input=input/example.csv --output=tmp/fine_results.csv --finegrain

Alternatively, you can run for coarse grain results using the following command:

    Rscript ./cps_run.R --input=input/example.csv --output=tmp/coarse_results.csv --coarsegrain

you can try SVR model for coarse grain results using the following command:

    Rscript ./cps_run.R --input=input/example.csv --output=tmp/coarse_svr_results.csv --coarsegrain --svr

### Using the functions

You can also use the CPS functions in your code. First, you have to load the source for '*Cell Proportion System*' code.

    source("cps_functions.R")

Load expression data from CSV file respecting the format already mentioned.

    expression_matrix <- read.csv("./input/example.csv",row.names = 1)

Run CPS function specifying if 'coarse' or 'fine' grained predictions are desired

    prediction_data <-do_CPS(expression_matrix, type="coarse", model="glmnet")

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

For the fine-grain case you just need to use:

    prediction_data <-do_CPS(expression_matrix, type="fine", model="glmnet")

and the function should return a data frame like this:

|                         |"Sample_1"         |"Sample_2"	    |
|-------------------------|-------------------|---------------------|
|"memory.B.cells"         |0.125474944057753  |0.0721011833308564   |
|"naive.B.cells"          |0.145859491671783  |0.138665106826305    |
|"memory.CD4.T.cells"     |0.196207716894681  |0.176457029951721    |
|"naive.CD4.T.cells"      |0.101451485325692  |0.115091362058503    |
|"regulatory.T.cells"     |0.0376165732796526 |0.0580151815019477   |
|"memory.CD8.T.cells"     |0.0820186245269355 |0.0548934376529131   |
|"naive.CD8.T.cells"      |0.0282401513096551 |-0.00500615193994581 |
|"NK.cells"               |0.0400782234299241 |0.0464635269670312   |
|"neutrophils"            |0.112686587369957  |0.116577437456931    |
|"monocytes"              |-0.0206236752459919|0.00714755867267107  |
|"myeloid.dendritic.cells"|0.0498252850971034 |0.0658293719429775   |
|"macrophages"            |0.20976358446132   |0.230940685578206    |
|"fibroblasts"            |0.122453816592088  |0.12388015785598     |
|"endothelial.cells"      |0.19206419268053   |0.183292376377571    |


The expected output for all the models is located in the `output/` directory

## Running time

Average run time for the `input/example.csv` dataset on a AMD Ryzen 7 with 32Gb is

svr    : 0m18,122s

glmnet : 0m2.,078s
