#' @export
dgcode <- function(data_list,
                   config = gcode::extract_config(verbose = F),
                   transfer = gcode::extract_transfer_framework(verbose = F),
                   recover = gcode::extract_recovery_framework(verbose = F),
                   join = gcode::extract_join_framework(verbose=F)
){
  
  copy_data_list <- data_list
  
  dynamic_list_of_models <- list()
  
  for (i in 1:config$bootstrap){
    
    print(paste("Bootstrap #    ",i,"   ...",sep=""))
    
    gcode.model <- gcode::gcode(data_list = data_list, config = config, join = join)
    dynamic_list_of_models <- c(dynamic_list_of_models,list(gcode.model))
    
    data_list <- lapply(c(1:length(data_list)),function(X){
      data_list[[X]] - t(gcode.model$main.parameters$alpha[[join$alpha[X]]])%*%gcode.model$main.code$code[[join$code[X]]]%*%t(gcode.model$main.parameters$beta[[join$beta[X]]])
    })
    
  }
  
  return(dynamic_list_of_models)
}
