#' Recover missing data points via imputation or prediction
#'
#' Main function to recover missing points. Used internally, but can be re-purposed by the user.
#'
#' @param data_list A list of datasets (matrix or tensor etc.)
#' @param code Code parameters from gcode
#' @param main.parameters Main parameters from gcode
#' @param config Configuration parameters from gcode
#' @param recover Recover list from gcode
#'
#' @return  Recovered data from imputation or prediction, with the design matrices and any user input parameters and functions
#' @export
recover_points <- function(data_list,
                           main.code,
                           main.parameters,
                           config,
                           join,
                           recover,
                           references,
                           iter){
  
  for (i in which(!is.null(recover$design.list))){
    
    missing_points <- which(as.matrix(recover$design.list[[i]]==1),arr.ind = T)
    row_with_missing_points <- unique(missing_points[,1])
    col_with_missing_points <- unique(missing_points[,2])
    
    Y <- data_list[[join$complete$data_list[i]]][,col_with_missing_points]
    Ytrain <- Y[-row_with_missing_points,]
    
    for (j in c(1:length(join$complete$data_list))){
      
      X <- (data_list[[join$complete$data_list[j]]][,-col_with_missing_points]%*%main.parameters$beta_sample[[join$complete$beta_sample[j]]][-col_with_missing_points,])
      Xtest <- cbind(1,X)[row_with_missing_points,]
      Xtrain <- cbind(1,X[-row_with_missing_points,])
      
      gcode.pred <- Xtest%*%MASS::ginv(t(Xtrain)%*%Xtrain)%*%t(Xtrain)%*%Ytrain
      glmnet.pred <- do.call('cbind',lapply(c(1:dim(Ytrain)[2]),function(X){
        glmnet.model <- glmnet::glmnet(x = scale(Xtrain), y = Ytrain[,X], family = "gaussian")
        rowMeans(predict(glmnet.model,newx = Xtest))
      }))
      knn.pred <- do.call('cbind',lapply(c(1:dim(Ytrain)[2]),function(X){
        FNN::knn.reg(train = Xtrain, test = Xtest, y = Ytrain[,X], k = 15)$pred
      }))
      rf.pred <- do.call('cbind',lapply(c(1:dim(Ytrain)[2]),function(X){
        ranger.model <- ranger::ranger(Ytrain[,X]~.,data=data.frame(x=Xtrain), num.trees = 15)
        predict(ranger.model,data.frame(x=Xtest))$predictions
      }))
      
      sample.pred <- (glmnet.pred+gcode.pred+knn.pred+rf.pred)/4
      data_list[[join$complete$data_list[i]]][row_with_missing_points,col_with_missing_points] <- sample.pred
      print(mean(abs(Ytrue-data_list[[join$complete$data_list[i]]][row_with_missing_points,col_with_missing_points])))
      
    }
    
  }
  
  return(list(
              data_list=data_list))
}
