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
                           references){
  
  for (i in which(do.call('c',lapply(recover$design.list,function(X){sum(X)>1})))){
    
    missing_points <- which(as.matrix(recover$design.list[[i]]==1),arr.ind = T)
    row_with_missing_points <- unique(missing_points[,1])
    column_with_missing_points <- unique(missing_points[,2])
    
    pred <- 0
    for (j in c(1:length(join$complete$data_list))[-i]){
      pred <- pred + (data_list[[join$complete$data_list[j]]] - main.parameters$intercept[[join$complete$data_list[j]]])%*%(main.parameters$beta_sample[[join$complete$beta_sample[j]]])%*%MASS::ginv(t((main.parameters$beta_sample[[join$complete$beta_sample[j]]]))%*%((main.parameters$beta_sample[[join$complete$beta_sample[j]]])))%*%t(main.parameters$beta_sample[[join$complete$beta_sample[i]]]) + main.parameters$intercept[[join$complete$data_list[i]]]
    }
    
    pred <- pred / c(length(join$complete$data_list) - 1 )
    
    data_list[[join$complete$data_list[i]]][row_with_missing_points,column_with_missing_points] <- pred[row_with_missing_points,column_with_missing_points]
    
  }
  
  return(list(
              data_list=data_list))
}
