#coarsegrained submission

library(glmnet)

# Read in the round and sub-Challenge-specific input file 
## listing each of the datasets
print(list.files())
print(getwd())

#input_df <- readr::read_csv("input/input.csv")

## Extract the names of each dataset

#dataset_names <- input_df$dataset.name

## Extract the names of the expression files that use 
## ENSEMBLE symbols as gene identifiers

expression_files <- "example.csv"

## Form the paths of the expression files
expression_paths <- paste0("input/", expression_files)


load("featureMaster.rdata")

#Function to make sample Normalization

Rhino_norm= function(l){
  
  #lr= rank(l*-1,ties.method="average")
  lr= rank(l*-1,ties.method="max")
  
  lr=log2(lr)
  lr= lr*-1
  lr= lr+ log2(length(lr))
  
  
  #lr[lr==min(lr)]=0
  res=lr
  return(res)
}

expression_matrix <- read.csv(expression_paths,row.names = 1)

if(any(!inALL %in% rownames(expression_matrix))){
  Fs=inALL[!inALL %in% rownames(expression_matrix)]
  Add_matrix= matrix(NA,nrow=length(Fs),ncol=ncol(expression_matrix))
  rownames(Add_matrix)=Fs; colnames(Add_matrix)= colnames(expression_matrix)
  expression_matrix=rbind(expression_matrix,Add_matrix)
  
}

expression_matrix= expression_matrix[inALL,]
expression_matrix=log2(expression_matrix+1)
expression_matrix[is.na(expression_matrix)]=log2(1)

X= apply(expression_matrix,2,Rhino_norm)
expression_matrix=X

#RUN CPS (Cell proportion system)

do_CPS <- function(expression_matrix, type= NA){
  
  if(type!= "coarse" & type !="fine"){stop("type should be 'coarse' or 'fine'")}
  
  
  #COURSE GRAIN
  if(type=="coarse"){  
    
    load("coarsegrained.rdata")
    p_B=predict(COURSEGRAIN$cvfit_B,newx=t(expression_matrix),s="lambda.1se")
    if(any(p_B<0)){ if(!all(p_B<=0))
    {p_B[p_B<0]=0}}
    
    p_CD4T=predict(COURSEGRAIN$cvfit_CD4T,newx=t(expression_matrix),s="lambda.1se")
    if(any(p_CD4T<0)){ if(!all(p_CD4T<=0))
    {p_CD4T[p_CD4T<0]=0}}
    
    p_CD8T=predict(COURSEGRAIN$cvfit_CD8T,newx=t(expression_matrix),s="lambda.1se")
    if(any(p_CD8T<0)){ if(!all(p_CD8T<=0))
    {p_CD8T[p_CD8T<0]=0}}
    
    p_NK=predict(COURSEGRAIN$cvfit_NK,newx=t(expression_matrix),s="lambda.1se")
    if(any(p_NK<0)){ if(!all(p_NK<=0))
    {p_NK[p_NK<0]=0}}
    
    p_Neutro=predict(COURSEGRAIN$cvfit_Neutro,newx=t(expression_matrix),s="lambda.1se")
    if(any(p_Neutro<0)){ if(!all(p_Neutro<=0))
    {p_Neutro[p_Neutro<0]=0}}
    
    p_Mono=predict(COURSEGRAIN$cvfit_Mono,newx=t(expression_matrix),s="lambda.1se")
    if(any(p_Mono<0)){ if(!all(p_Mono<=0))
    {p_Mono[p_Mono<0]=0}}
    
    p_Fibro=predict(COURSEGRAIN$cvfit_Fibro,newx=t(expression_matrix),s="lambda.1se")
    if(any(p_Fibro<0)){ if(!all(p_Fibro<=0))
    {p_Fibro[p_Fibro<0]=0}}
    
    p_Endo=predict(COURSEGRAIN$cvfit_Endo,newx=t(expression_matrix),s="lambda.1se")
    if(any(p_Endo<0)){ if(!all(p_Endo<=0))
    {p_B[p_Endo<0]=0}}
    
    result_matrix= data.frame(B.cells=as.numeric(p_B),
                              CD4.T.cells=as.numeric(p_CD4T),
                              CD8.T.cells=as.numeric(p_CD8T),
                              NK.cells=as.numeric(p_NK),
                              neutrophils=as.numeric(p_Neutro),
                              monocytic.lineage=as.numeric(p_Mono),
                              fibroblasts=as.numeric(p_Fibro),
                              endothelial.cells=as.numeric(p_Endo))
    rownames(result_matrix)= colnames(expression_matrix)
    result_matrix=as.data.frame(t(result_matrix))
    
    
    ## Create the directory the output will go into
    dir.create("output")
    ## Write result into output directory
    write.csv(result_matrix, "output/predictions_Coarsegrain.csv",row.names=TRUE)
  }
  
  if(type=="fine"){
    
    load("finegrained.rdata")
    p_MBC=predict(FINEGRAIN$cvfit_MBC,newx=t(expression_matrix),s="lambda.1se")
    p_NBC=predict(FINEGRAIN$cvfit_NBC,newx=t(expression_matrix),s="lambda.1se")
    p_MCD4T=predict(FINEGRAIN$cvfit_MCD4T,newx=t(expression_matrix),s="lambda.1se")
    p_NCD4T=predict(FINEGRAIN$cvfit_NCD4T,newx=t(expression_matrix),s="lambda.1se")
    p_Treg=predict(FINEGRAIN$cvfit_Treg,newx=t(expression_matrix),s="lambda.1se")
    p_MCD8=predict(FINEGRAIN$cvfit_MCD8,newx=t(expression_matrix),s="lambda.1se")
    p_NCD8=predict(FINEGRAIN$cvfit_NCD8,newx=t(expression_matrix),s="lambda.1se")
    p_MDC=predict(FINEGRAIN$cvfit_MDC,newx=t(expression_matrix),s="lambda.1se")
    p_Macro=predict(FINEGRAIN$cvfit_Macro,newx=t(expression_matrix),s="lambda.1se")
    p_Fibro=predict(FINEGRAIN$cvfit_Fibro,newx=t(expression_matrix),s="lambda.1se")
    p_Endo=predict(FINEGRAIN$cvfit_Endo,newx=t(expression_matrix),s="lambda.1se")
    p_Neutro=predict(FINEGRAIN$cvfit_Neutro,newx=t(expression_matrix),s="lambda.1se")
    p_NK=predict(FINEGRAIN$cvfit_NK,newx=t(expression_matrix),s="lambda.1se")
    p_Monocyte=predict(FINEGRAIN$cvfit_Monocyte,newx=t(expression_matrix),s="lambda.1se")
    
    result_matrix= data.frame(
      memory.B.cells=as.numeric(p_MBC),
      naive.B.cells= as.numeric(p_NBC),
      memory.CD4.T.cells= as.numeric(p_MCD4T),
      naive.CD4.T.cells=as.numeric(p_NCD4T),
      regulatory.T.cells= as.numeric(p_Treg),
      memory.CD8.T.cells= as.numeric(p_MCD8),
      naive.CD8.T.cells= as.numeric(p_NCD8),
      NK.cells= as.numeric(p_NK),
      neutrophils=as.numeric(p_Neutro),
      monocytes=as.numeric(p_Monocyte),
      myeloid.dendritic.cells=as.numeric(p_MDC),
      macrophages=as.numeric(p_Macro),
      fibroblasts=as.numeric(p_Fibro),
      endothelial.cells=as.numeric(p_Endo))
    
    rownames(result_matrix)= colnames(expression_matrix)
    result_matrix=as.data.frame(t(result_matrix))
    
    ## Create the directory the output will go into
    dir.create("output")
    ## Write result into output directory
    write.csv(result_matrix, "output/predictions_Finegrain.csv",row.names=TRUE)
    
  }  
}

do_CPS(expression_matrix,type="coarse")
