#' @export
initialise.gcode <- function(data_list,
                             config,
                             transfer,
                             join
){
  
  if (config$verbose){
    print(c(paste("Initialising alpha with : ",config$init[1],sep=""),
            paste("Initialising beta with : ",config$init[2],sep=""))
    )
  }
  
  data_list <- lapply(data_list,function(X){
    X <- gcode::prepare_data(X,
                        log=config$transform$log,
                        center=config$transform$center,
                        scale.norm=config$transform$norm
                        )
  })
  
  main.code <- list(code=lapply(unique(join$complete$code),function(X){NULL}),encode=lapply(unique(join$complete$code),function(X){NULL}),incode=lapply(unique(join$complete$incode),function(X){NULL}))
  main.parameters <- list(alpha = lapply(unique(join$complete$alpha),function(X){NULL}), beta = lapply(unique(join$complete$beta),function(X){NULL}),
                          alpha.code = lapply(unique(join$complete$alpha.code),function(X){NULL}), beta.code = lapply(unique(join$complete$beta.code),function(X){NULL}),
                          intercept = lapply(unique(1:length(join$complete$data_list)),function(X){NULL}))
  
  config.code <- config
  config.code$i_dim <- config.code$k_dim
  config.code$j_dim <- config.code$k_dim
  
  for (i in 1:length(join$complete$data_list)){
    
    if (is.null(main.parameters$alpha[[join$complete$alpha[i]]])){
      if (!is.null(transfer$main.parameters$alpha[[join$complete$alpha[i]]])){
        main.parameters$alpha[[join$complete$alpha[i]]] <- transfer$main.parameters$alpha[[join$complete$alpha[i]]]
      } else {
        main.parameters$alpha[[join$complete$alpha[i]]] <- as.matrix(initialise.parameters(x = data_list[[join$complete$data_list[i]]],config = config, param.type = "alpha"))
      }
    }
    
    if (is.null(main.parameters$beta[[join$complete$beta[i]]])){
      if (!is.null(transfer$main.parameters$beta[[join$complete$beta[i]]])){
        main.parameters$beta[[join$complete$beta[i]]] <- transfer$main.parameters$beta[[join$complete$beta[i]]]
      } else {
        main.parameters$beta[[join$complete$beta[i]]] <- as.matrix(initialise.parameters(x = data_list[[join$complete$data_list[i]]],config = config, param.type = "beta"))
      }
    }
    
    
    
    
    if (is.null(main.parameters$alpha.code[[join$complete$alpha.code[i]]])){
      if (!is.null(transfer$main.parameters$alpha.code[[join$complete$alpha.code[i]]])){
        main.parameters$alpha.code[[join$complete$alpha.code[i]]] <- transfer$main.parameters$alpha.code[[join$complete$alpha.code[i]]]
      } else {
        main.parameters$alpha.code[[join$complete$alpha.code[i]]] <- t(as.matrix(initialise.parameters(x = diag(config$i_dim),config = config.code, param.type = "alpha")))
      }
    }
    
    if (is.null(main.parameters$beta.code[[join$complete$beta.code[i]]])){
      if (!is.null(transfer$main.parameters$beta.code[[join$complete$beta.code[i]]])){
        main.parameters$beta.code[[join$complete$beta.code[i]]] <- transfer$main.parameters$beta.code[[join$complete$beta.code[i]]]
      } else {
        main.parameters$beta.code[[join$complete$beta.code[i]]] <- t(as.matrix(initialise.parameters(x = diag(config$j_dim),config = config.code, param.type = "beta")))
      }
    }
    
    
    
    
    
    if (is.null(main.code$encode[[join$complete$code[i]]])){
      if (!is.null(transfer$main.code$code[[join$complete$code[i]]])){
        main.code$encode[[join$complete$code[i]]] <- transfer$main.code$encode[[join$complete$code[i]]]
      } else {
        main.code$encode[[join$complete$code[i]]] <- main.parameters$alpha[[join$complete$alpha[i]]]%*%as.matrix(data_list[[join$complete$data_list[i]]])%*%main.parameters$beta[[join$complete$beta[i]]]
      }
    }

    if (is.null(main.code$code[[join$complete$code[i]]])){
      if (!is.null(transfer$main.code$code[[join$complete$code[i]]])){
        main.code$code[[join$complete$code[i]]] <- transfer$main.code$code[[join$complete$code[i]]]
      } else {
        main.code$code[[join$complete$code[i]]] <- as.matrix(MASS::ginv(main.parameters$alpha[[join$complete$alpha[i]]]%*%t(main.parameters$alpha[[join$complete$alpha[i]]]))%*%main.code$encode[[join$complete$code[i]]]%*%MASS::ginv(t(main.parameters$beta[[join$complete$beta[i]]])%*%main.parameters$beta[[join$complete$beta[i]]]))
      }
    }
    
    if (is.null(main.code$incode[[join$complete$incode[i]]])){
      if (!is.null(transfer$main.code$code[[join$complete$code[i]]])){
        main.code$incode[[join$complete$incode[i]]] <- transfer$main.code$incode[[join$complete$incode[i]]]
      } else {
        main.code$incode[[join$complete$incode[i]]] <- (t(main.parameters$alpha.code[[join$complete$alpha.code[i]]])%*%as.matrix(main.code$code[[join$complete$code[i]]])%*%t(main.parameters$beta.code[[join$complete$beta.code[i]]]))
      }
    }
    
  
    
    if (is.null(main.parameters$intercept[[join$complete$data_list[i]]])){
      if (!is.null(transfer$main.parameters$intercept[[join$complete$data_list[i]]])){
        main.parameters$intercept[[join$complete$data_list[i]]] <- transfer$main.parameters$intercept[[join$complete$data_list[i]]]
      } else {
        main.parameters$intercept[[join$complete$data_list[i]]] <- colMeans(data_list[[join$complete$data_list[i]]] - t(main.parameters$alpha[[join$complete$alpha[i]]])%*%main.code$code[[join$complete$code[i]]]%*%t(main.parameters$beta[[join$complete$beta[i]]]))
      }
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
      (irlba::irlba(as.matrix(x), nv = config$j_dim, maxit = 15)$v)
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
      t(irlba::irlba(as.matrix(x), nu = config$i_dim, maxit = 15)$u)
    } else if (config$init[[1]]=="rsvd") {
      t(rsvd::rsvd(as.matrix(x), nu = config$i_dim)$u)
    } 
  }
  
  return(main.param)
  
}