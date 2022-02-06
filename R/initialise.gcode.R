

#' @export
initialise.gcode <- function(data_list,
                             config,
                             transfer,
                             join
){
  
  if (config$verbose){
    print(paste("Initialising data with : ",config$init,sep=""))
  }
  
  main.code <- list(code=lapply(unique(join$code),function(X){NULL}),encode=lapply(unique(join$code),function(X){NULL}))
  main.parameters <- list(alpha = lapply(unique(join$alpha),function(X){NULL}), beta = lapply(unique(join$beta),function(X){NULL}))
  
  for (iter in 1:2){
    for (i in 1:length(data_list)){

      if (is.null(main.parameters$alpha[[join$alpha[i]]])){
        if (!is.null(transfer$main.parameters$alpha[[join$alpha[i]]])){
          main.parameters$alpha[[join$alpha[i]]] <- transfer$main.parameters$alpha[[join$alpha[i]]]
        } else {
          main.parameters$alpha[[join$alpha[i]]] <- as.matrix(initialise.parameters(x = data_list[[i]],config = config, param.type = "alpha"))
        }
      }
      if (is.null(main.parameters$beta[[join$beta[i]]])){
        if (!is.null(transfer$main.parameters$beta[[join$beta[i]]])){
          main.parameters$beta[[join$beta[i]]] <- transfer$main.parameters$beta[[join$beta[i]]]
        } else {
          main.parameters$beta[[join$beta[i]]] <- as.matrix(initialise.parameters(x = data_list[[i]],config = config, param.type = "beta"))
        }
        
      }
      
      internal.parameters <- list(alpha=main.parameters$alpha[[join$alpha[i]]],
                                  beta=main.parameters$beta[[join$beta[i]]])
      
      if (is.null(main.code$code[[join$code[i]]]) | is.null(main.code$encode[[join$code[i]]])){
        if (!is.null(transfer$main.code$code[[join$code[i]]])){
          main.code$encode[[join$code[i]]] <- transfer$main.code$encode[[join$code[i]]]
          main.code$code[[join$code[i]]] <- transfer$main.code$code[[join$code[i]]]
        } else {
          main.code$encode[[join$code[i]]] <- main.parameters$alpha[[join$alpha[i]]]%*%as.matrix(data_list[[i]])%*%main.parameters$beta[[join$beta[i]]]
          main.code$code[[join$code[i]]] <- MASS::ginv(main.parameters$alpha[[join$alpha[i]]]%*%t(main.parameters$alpha[[join$alpha[i]]]))%*%main.code$encode[[join$code[i]]]%*%MASS::ginv(t(main.parameters$beta[[join$beta[i]]])%*%main.parameters$beta[[join$beta[i]]])
        }
      }
      
      internal.code <- list(encode=main.code$encode[[join$code[i]]],
                            code=main.code$code[[join$code[i]]])
      
      return_update <- update_set(x = as.matrix(data_list[[i]]),
                                  main.parameters = internal.parameters,
                                  main.code = internal.code, 
                                  config = config,
                                  fix = transfer$fix
      )
      
      main.parameters$alpha[[join$alpha[i]]] <- return_update$main.parameters$alpha
      main.parameters$beta[[join$beta[i]]] <- return_update$main.parameters$beta

      main.code$code[[join$code[i]]] <- return_update$main.code$code
      main.code$encode[[join$code[i]]] <- return_update$main.code$encode
      
    }
  }
  
  return(
    list(
      main.parameters = main.parameters,
      main.code = main.code
    )
  )
  
}


#' @export
initialise.parameters <- function(x,config,param.type){
  
  if (param.type=="beta"){
    main.param <- if (config$init[[2]]=="rnorm"){
      array(rnorm(dim(x)[2]*config$j_dim),dim=c(dim(x)[2],config$j_dim))
    } else if (config$init[[2]]=="runif"){
      array(runif(dim(x)[2]*config$j_dim),dim=c(dim(x)[2],config$j_dim))
    } else if (config$init[[2]]=="irlba") {
      (irlba::irlba(as.matrix(x), nv = config$j_dim, maxit = 50)$v)
    } else if (config$init[[2]]=="rsvd") {
      (rsvd::rsvd(as.matrix(x), nv = config$j_dim)$v)
    } 
  }
  
  if (param.type=="alpha"){
    main.param <- if (config$init[[1]]=="rnorm") {
      array(rnorm(config$i_dim*dim(x)[1]),dim=c(config$i_dim,dim(x)[1]))
    } else if (config$init[[1]]=="runif") {
      array(runif(config$i_dim*dim(x)[1]),dim=c(config$i_dim,dim(x)[1]))
    } else if (config$init[[1]]=="irlba") {
      t(irlba::irlba(as.matrix(x), nu = config$i_dim, maxit = 50)$u)
    } else if (config$init[[1]]=="rsvd") {
      t(rsvd::rsvd(as.matrix(x), nu = config$i_dim)$u)
    } 
  }
  
  return(main.param)
  
}