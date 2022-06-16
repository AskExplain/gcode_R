#' @export
dgcode <- function(data_list,
                   config = gcode::extract_config(verbose = F),
                   transfer = gcode::extract_transfer_framework(verbose = F),
                   recover = gcode::extract_recovery_framework(verbose = F),
                   join = gcode::extract_join_framework(verbose=F),
                   references = gcode::extract_references_framework(F)
){
  
  combin <- cbind(rep(config$i_dim,length(config$j_dim)),rep(config$j_dim,length(config$i_dim)))
  
  main_data_list <- data_list
  main_join <- join
  
  config$i_dim <- combin[1,1]
  config$j_dim <- combin[1,2]
  config$verbose <- F
  
  dynamic_list_of_models <- list()
  
  gcode.model <- gcode::gcode(data_list = data_list, config = config, join = join, transfer = transfer, recover = recover, references = references)
  dynamic_list_of_models <- c(dynamic_list_of_models,list(gcode.model))
  
  decoded <- c(lapply(c(1:length(main_join$complete$data_list)),function(X){
    t(gcode.model$main.parameters$alpha_sample[[main_join$complete$alpha_sample[X]]])%*%gcode.model$main.code$code[[main_join$complete$code[X]]]%*%t(gcode.model$main.parameters$beta_sample[[main_join$complete$beta_sample[X]]]) + gcode.model$main.parameters$intercept[[main_join$complete$data_list[X]]]
  }))
  
  data_list <- c(decoded,main_data_list)
  
  join$complete <- lapply(main_join$complete,function(X){c(X,X+max(X))})
  join$complete$alpha_signal[1:length(join$complete$alpha_signal)] <- 1

  for (i in 2:dim(combin)[1]){
    
    config$i_dim <- combin[i,1]
    config$j_dim <- combin[i,2]
    
    config$verbose <- F
    print(paste("Bootstrap: iteration   ",i,"   of   ",dim(combin)[1],"   ...",sep=""))
    
    gcode.model <- gcode::gcode(data_list = data_list, config = config, join = join, transfer = transfer, recover = recover, references = references)
    dynamic_list_of_models <- c(dynamic_list_of_models,list(gcode.model))
    
    decoded <- c(lapply(c(1:length(main_join$complete$data_list)),function(X){
      t(gcode.model$main.parameters$alpha_sample[[main_join$complete$alpha_sample[X]]])%*%gcode.model$main.code$code[[main_join$complete$code[X]]]%*%t(gcode.model$main.parameters$beta_sample[[main_join$complete$beta_sample[X]]]) + gcode.model$main.parameters$intercept[[main_join$complete$data_list[X]]]
    }))
    
    data_list <- c(decoded,main_data_list)
    
    join$complete <- lapply(main_join$complete,function(X){c(X,X+max(X))})
    join$complete$alpha_signal[1:length(join$complete$alpha_signal)] <- 1

  }
  
  return(dynamic_list_of_models)
}
