# Dream Deconvolution Challenge Team DA_505
Team DA_505 submissions to DREAM Challenge Tumor Deconvolution

### Installation
1. Clone the GitHub repo.
```
git clone https://github.com/martinguerrero89/Dream_Deconv_Challenge_Team_DA505.git
```   
2. Install the following requirements.
```
R >= 4.0.0
glmnet
```
3. In R Enter the workspace.
```
setwd("./Dream_Deconv_Challenge_Team_DA505")
```
### Usage

Save your input expression matrix in the 'input' folder and open it in R.
The expression data should be a comma separated file with columns associated to sample ID and rows to genes specified using ENSG codes.

```
expression_matrix <- read.csv("./input/example.csv",row.names = 1)
```

Load the source for 'Cell Proportion System' code

```
source("RUN_CPS.R")
```

Run CPS function specifying if 'coarse' or 'fine' grained predictions are desired

```
do_CPS(expression_matrix, type="coarse")
```

Predictions will be in a 'output' folder

