
#' Title
#'
#' @param l 
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param expresion_matrix 
#'
#' @return
#' @export
#'
#' @examples
preprocess_em<- function(expression_matrix){
  load("models/featureMaster.rdata")
  if(any(!inALL %in% rownames(expression_matrix))){
    Fs=inALL[!inALL %in% rownames(expression_matrix)]
    Add_matrix= matrix(NA,nrow=length(Fs),ncol=ncol(expression_matrix))
    rownames(Add_matrix)=Fs; colnames(Add_matrix)= colnames(expression_matrix)
    expression_matrix=rbind(expression_matrix,Add_matrix)
  }
  expression_matrix <- expression_matrix[inALL,]
  expression_matrix <- log2(expression_matrix+1)
  expression_matrix[is.na(expression_matrix)]=log2(1)
  X <- apply(expression_matrix,2,Rhino_norm)
  expression_matrix<-X
}

#' Title
#'
#' @param expression_matrix 
#'
#' @return
#' @export
#'
#' @examples
do_CPS_coarse_glmnet <- function(expression_matrix){
    expression_matrix <- preprocess_em(expression_matrix)
    load("models/glmnet_model_coarsegrained.rdata")
    p_B = predict(COURSEGRAIN$cvfit_B,
                  newx = t(expression_matrix),
                  s = "lambda.1se")
    if (any(p_B < 0)) {
      if (!all(p_B <= 0))
      {
        p_B[p_B < 0] = 0
      }
    }
    
    p_CD4T = predict(COURSEGRAIN$cvfit_CD4T,
                     newx = t(expression_matrix),
                     s = "lambda.1se")
    if (any(p_CD4T < 0)) {
      if (!all(p_CD4T <= 0))
      {
        p_CD4T[p_CD4T < 0] = 0
      }
    }
    
    p_CD8T = predict(COURSEGRAIN$cvfit_CD8T,
                     newx = t(expression_matrix),
                     s = "lambda.1se")
    if (any(p_CD8T < 0)) {
      if (!all(p_CD8T <= 0))
      {
        p_CD8T[p_CD8T < 0] = 0
      }
    }
    
    p_NK = predict(COURSEGRAIN$cvfit_NK,
                   newx = t(expression_matrix),
                   s = "lambda.1se")
    if (any(p_NK < 0)) {
      if (!all(p_NK <= 0))
      {
        p_NK[p_NK < 0] = 0
      }
    }
    
    p_Neutro = predict(COURSEGRAIN$cvfit_Neutro,
                       newx = t(expression_matrix),
                       s = "lambda.1se")
    if (any(p_Neutro < 0)) {
      if (!all(p_Neutro <= 0))
      {
        p_Neutro[p_Neutro < 0] = 0
      }
    }
    
    p_Mono = predict(COURSEGRAIN$cvfit_Mono,
                     newx = t(expression_matrix),
                     s = "lambda.1se")
    if (any(p_Mono < 0)) {
      if (!all(p_Mono <= 0))
      {
        p_Mono[p_Mono < 0] = 0
      }
    }
    
    p_Fibro = predict(COURSEGRAIN$cvfit_Fibro,
                      newx = t(expression_matrix),
                      s = "lambda.1se")
    if (any(p_Fibro < 0)) {
      if (!all(p_Fibro <= 0))
      {
        p_Fibro[p_Fibro < 0] = 0
      }
    }
    
    p_Endo = predict(COURSEGRAIN$cvfit_Endo,
                     newx = t(expression_matrix),
                     s = "lambda.1se")
    if (any(p_Endo < 0)) {
      if (!all(p_Endo <= 0))
      {
        p_B[p_Endo < 0] = 0
      }
    }
    
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
}  
  
#' Title
#'
#' @param expression_matrix 
#'
#' @return
#' @export
#'
#' @examples
do_CPS_coarse_svr <- function(expression_matrix) {
  expression_matrix <- preprocess_em(expression_matrix)
  load("models/svr_model_coarsegrained_B.cells.rdata")
  best_features <-
    read.csv("models/svr_model_coarsegrained_B.cells_selected_variables",
             header = FALSE)
  p_B = predict(COURSEGRAIN_SVR,
                newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),]))
  if (any(p_B < 0)) {
    if (!all(p_B <= 0))
    {
      p_B[p_B < 0] = 0
    }
  }
  load("models/svr_model_coarsegrained_CD4.T.cells.rdata")
  best_features <-
    read.csv("models/svr_model_coarsegrained_CD4.T.cells_selected_variables",
             header = FALSE)
  
  p_CD4T = predict(COURSEGRAIN_SVR,
                   newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),]))
  if (any(p_CD4T < 0)) {
    if (!all(p_CD4T <= 0))
    {
      p_CD4T[p_CD4T < 0] = 0
    }
  }
  load("models/svr_model_coarsegrained_CD8.T.cells.rdata")
  best_features <-
    read.csv("models/svr_model_coarsegrained_CD8.T.cells_selected_variables",
             header = FALSE)
  p_CD8T = predict(COURSEGRAIN_SVR,
                   newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),]))
  if (any(p_CD8T < 0)) {
    if (!all(p_CD8T <= 0))
    {
      p_CD8T[p_CD8T < 0] = 0
    }
  }
  load("models/svr_model_coarsegrained_NK.cells.rdata")
  best_features <-
    read.csv("models/svr_model_coarsegrained_NK.cells_selected_variables",
             header = FALSE)
  p_NK = predict(COURSEGRAIN_SVR,
                 newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),]))
  if (any(p_NK < 0)) {
    if (!all(p_NK <= 0))
    {
      p_NK[p_NK < 0] = 0
    }
  }
  load("models/svr_model_coarsegrained_neutrophils.rdata")
  best_features <-
    read.csv("models/svr_model_coarsegrained_neutrophils_selected_variables",
             header = FALSE)
  p_Neutro = predict(COURSEGRAIN_SVR,
                     newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),]))
  if (any(p_Neutro < 0)) {
    if (!all(p_Neutro <= 0))
    {
      p_Neutro[p_Neutro < 0] = 0
    }
  }
  load("models/svr_model_coarsegrained_monocytic.lineage.rdata")
  best_features <-
    read.csv(
      "models/svr_model_coarsegrained_monocytic.lineage_selected_variables",
      header = FALSE
    )
  p_Mono = predict(COURSEGRAIN_SVR,
                   newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),]))
  if (any(p_Mono < 0)) {
    if (!all(p_Mono <= 0))
    {
      p_Mono[p_Mono < 0] = 0
    }
  }
  load("models/svr_model_coarsegrained_fibroblasts.rdata")
  best_features <-
    read.csv("models/svr_model_coarsegrained_fibroblasts_selected_variables",
             header = FALSE)
  p_Fibro = predict(COURSEGRAIN_SVR,
                    newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),]))
  if (any(p_Fibro < 0)) {
    if (!all(p_Fibro <= 0))
    {
      p_Fibro[p_Fibro < 0] = 0
    }
  }
  
  load("models/svr_model_coarsegrained_endothelial.cells.rdata")
  best_features <-
    read.csv(
      "models/svr_model_coarsegrained_endothelial.cells_selected_variables",
      header = FALSE
    )
  p_Endo = predict(COURSEGRAIN_SVR,
                   newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),]))
  if (any(p_Endo < 0)) {
    if (!all(p_Endo <= 0))
    {
      p_B[p_Endo < 0] = 0
    }
  }
  
  result_matrix = data.frame(
    B.cells = as.numeric(p_B),
    CD4.T.cells = as.numeric(p_CD4T),
    CD8.T.cells = as.numeric(p_CD8T),
    NK.cells = as.numeric(p_NK),
    neutrophils = as.numeric(p_Neutro),
    monocytic.lineage = as.numeric(p_Mono),
    fibroblasts = as.numeric(p_Fibro),
    endothelial.cells = as.numeric(p_Endo)
  )
  rownames(result_matrix) = colnames(expression_matrix)
  result_matrix = as.data.frame(t(result_matrix))
}  

#' Title
#'
#' @param expression_matrix 
#'
#' @return
#' @export
#'
#' @examples
do_CPS_fine_svr <- function(expression_matrix) {
  expression_matrix <- preprocess_em(expression_matrix)
  load("models/svr_model_finegrained_memory.B.cells.rdata")
  best_features = read.csv("models/svr_model_finegrained_memory.B.cells_selected_variables",header = FALSE)
  p_MBC = predict(FINEGRAIN_SVR,
                 newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
                 )
  load("models/svr_model_finegrained_naive.B.cells.rdata")
  best_features = read.csv("models/svr_model_finegrained_naive.B.cells_selected_variables",header = FALSE)
  p_NBC = predict(FINEGRAIN_SVR,
                  newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
                  )
  load("models/svr_model_finegrained_memory.CD4.T.cells.rdata")
  best_features = read.csv("models/svr_model_finegrained_memory.CD4.T.cells_selected_variables",header = FALSE)
  p_MCD4T = predict(FINEGRAIN_SVR,
                    newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
                    )
  load("models/svr_model_finegrained_naive.CD4.T.cells.rdata")
  best_features = read.csv("models/svr_model_finegrained_naive.CD4.T.cells_selected_variables",header = FALSE)
  p_NCD4T = predict(FINEGRAIN_SVR,
                    newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
                    )
  load("models/svr_model_finegrained_regulatory.T.cells.rdata")
  best_features = read.csv("models/svr_model_finegrained_regulatory.T.cells_selected_variables",header = FALSE)
  p_Treg = predict(FINEGRAIN_SVR,
                   newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
                    )
  load("models/svr_model_finegrained_memory.CD8.T.cells.rdata")
  best_features = read.csv("models/svr_model_finegrained_memory.CD8.T.cells_selected_variables",header = FALSE)
  p_MCD8 = predict(FINEGRAIN_SVR,
                   newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
  )
  load("models/svr_model_finegrained_naive.CD8.T.cells.rdata")
  best_features = read.csv("models/svr_model_finegrained_naive.CD8.T.cells_selected_variables",header = FALSE)
  p_NCD8 = predict(FINEGRAIN_SVR,
                   newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
  )
  load("models/svr_model_finegrained_myeloid.dendritic.cells.rdata")
  best_features = read.csv("models/svr_model_finegrained_myeloid.dendritic.cells_selected_variables",header = FALSE)
  p_MDC = predict(FINEGRAIN_SVR,
                  newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
  )
  load("models/svr_model_finegrained_macrophages.rdata")
  best_features = read.csv("models/svr_model_finegrained_macrophages_selected_variables",header = FALSE)
  p_Macro = predict(FINEGRAIN_SVR,
                    newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
  )
  load("models/svr_model_finegrained_fibroblasts.rdata")
  best_features = read.csv("models/svr_model_finegrained_fibroblasts_selected_variables",header = FALSE)
  p_Fibro = predict(FINEGRAIN_SVR,
                     newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
  )
  load("models/svr_model_finegrained_endothelial.cells.rdata")
  best_features = read.csv("models/svr_model_finegrained_endothelial.cells_selected_variables",header = FALSE)
  p_Endo = predict(FINEGRAIN_SVR,
                   newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
  )
  
  load("models/svr_model_finegrained_neutrophils.rdata")
  best_features = read.csv("models/svr_model_finegrained_neutrophils_selected_variables",header = FALSE)
  p_Neutro = predict(FINEGRAIN_SVR,
                     newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
  )
  
  load("models/svr_model_finegrained_NK.cells.rdata")
  best_features = read.csv("models/svr_model_finegrained_NK.cells_selected_variables",header = FALSE)
  p_NK = predict(FINEGRAIN_SVR,
                 newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
  )
  
  load("models/svr_model_finegrained_monocytes.rdata")
  best_features = read.csv("models/svr_model_finegrained_monocytes_selected_variables",header = FALSE)
  p_Monocyte = predict(FINEGRAIN_SVR,
                       newdata = t(expression_matrix[which(rownames(expression_matrix) %in% as.matrix(best_features)),])
  )
  
  result_matrix = data.frame(
    memory.B.cells = as.numeric(p_MBC),
    naive.B.cells = as.numeric(p_NBC),
    memory.CD4.T.cells = as.numeric(p_MCD4T),
    naive.CD4.T.cells = as.numeric(p_NCD4T),
    regulatory.T.cells = as.numeric(p_Treg),
    memory.CD8.T.cells = as.numeric(p_MCD8),
    naive.CD8.T.cells = as.numeric(p_NCD8),
    NK.cells = as.numeric(p_NK),
    neutrophils = as.numeric(p_Neutro),
    monocytes = as.numeric(p_Monocyte),
    myeloid.dendritic.cells = as.numeric(p_MDC),
    macrophages = as.numeric(p_Macro),
    fibroblasts = as.numeric(p_Fibro),
    endothelial.cells = as.numeric(p_Endo)
  )
  
  rownames(result_matrix) = colnames(expression_matrix)
  result_matrix = as.data.frame(t(result_matrix))
}


#' Title
#'
#' @param expression_matrix 
#'
#' @return
#' @export
#'
#' @examples
do_CPS_fine_glmnet <- function(expression_matrix) {
  expression_matrix <- preprocess_em(expression_matrix)
  load("models/glmnet_model_finegrained.rdata")
  p_MBC = predict(FINEGRAIN$cvfit_MBC,
                  newx = t(expression_matrix),
                  s = "lambda.1se")
  p_NBC = predict(FINEGRAIN$cvfit_NBC,
                  newx = t(expression_matrix),
                  s = "lambda.1se")
  p_MCD4T = predict(FINEGRAIN$cvfit_MCD4T,
                    newx = t(expression_matrix),
                    s = "lambda.1se")
  p_NCD4T = predict(FINEGRAIN$cvfit_NCD4T,
                    newx = t(expression_matrix),
                    s = "lambda.1se")
  p_Treg = predict(FINEGRAIN$cvfit_Treg,
                   newx = t(expression_matrix),
                   s = "lambda.1se")
  p_MCD8 = predict(FINEGRAIN$cvfit_MCD8,
                   newx = t(expression_matrix),
                   s = "lambda.1se")
  p_NCD8 = predict(FINEGRAIN$cvfit_NCD8,
                   newx = t(expression_matrix),
                   s = "lambda.1se")
  p_MDC = predict(FINEGRAIN$cvfit_MDC,
                  newx = t(expression_matrix),
                  s = "lambda.1se")
  p_Macro = predict(FINEGRAIN$cvfit_Macro,
                    newx = t(expression_matrix),
                    s = "lambda.1se")
  p_Fibro = predict(FINEGRAIN$cvfit_Fibro,
                    newx = t(expression_matrix),
                    s = "lambda.1se")
  p_Endo = predict(FINEGRAIN$cvfit_Endo,
                   newx = t(expression_matrix),
                   s = "lambda.1se")
  p_Neutro = predict(FINEGRAIN$cvfit_Neutro,
                     newx = t(expression_matrix),
                     s = "lambda.1se")
  p_NK = predict(FINEGRAIN$cvfit_NK,
                 newx = t(expression_matrix),
                 s = "lambda.1se")
  p_Monocyte = predict(FINEGRAIN$cvfit_Monocyte,
                       newx = t(expression_matrix),
                       s = "lambda.1se")
  
  result_matrix = data.frame(
    memory.B.cells = as.numeric(p_MBC),
    naive.B.cells = as.numeric(p_NBC),
    memory.CD4.T.cells = as.numeric(p_MCD4T),
    naive.CD4.T.cells = as.numeric(p_NCD4T),
    regulatory.T.cells = as.numeric(p_Treg),
    memory.CD8.T.cells = as.numeric(p_MCD8),
    naive.CD8.T.cells = as.numeric(p_NCD8),
    NK.cells = as.numeric(p_NK),
    neutrophils = as.numeric(p_Neutro),
    monocytes = as.numeric(p_Monocyte),
    myeloid.dendritic.cells = as.numeric(p_MDC),
    macrophages = as.numeric(p_Macro),
    fibroblasts = as.numeric(p_Fibro),
    endothelial.cells = as.numeric(p_Endo)
  )
  
  rownames(result_matrix) = colnames(expression_matrix)
  result_matrix = as.data.frame(t(result_matrix))
}

#' Title
#'
#' @param expression_matrix 
#' @param type 
#' @param model
#'
#' @return predictions
#' @export
#'
#' @examples
do_CPS <- function(expression_matrix, type= NA, model="glmnet"){
  if(type!= "coarse" & type !="fine"){stop("type should be 'coarse' or 'fine'")}
  if(model!= "svr" & type !="glmnet"){stop("model should be 'svr' or 'glmnet'")}
  
  if (type == "coarse") {
    if (model == "glmnet") predictions <- do_CPS_coarse_glmnet(expression_matrix)
    else predictions <- do_CPS_coarse_svr(expression_matrix)
    
  } else{
    if ( model == "glmnet") predictions <- do_CPS_fine_glmnet(expression_matrix)
    else predictions <- do_CPS_fine_svr(expression_matrix)
  }
  predictions
}  