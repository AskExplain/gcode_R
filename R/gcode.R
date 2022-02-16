#' Generative Encoding via Generalised Canonical Procrustes
#'
#' A method that uses a likelihood model to align multiple datasets via an encoding in a lower dimensional space. The parameters can be used to reduce either the feature or the sample dimensions into a smaller subspace for further embedding or prediction. To run as default, only a data list is required - please review the config parameters at gcode::extract_config(T)  .
#'
#' @param data_list List of data matrices of varying dimensionality. Attempts to find similarities among all datasets with a core structure.
#' @param config Configuration parameters (required, default provided)
#' @param transfer Transferring pre-trained model parameters (not required)
#' @param recover Important information used for prediction or imputation (not required)
#' @param join Join parameters that share the same axes to be similar (not required)
#'
#' @return Main parameters contains the learned model parameters. The alpha and beta matrix multiply example datasets x and y by, (K)(Y)(v) and (L)(X)(u). By multiplying with the parameter, the dimension of the samples and features can be dimensionally reduced for further visualisation analysis such as embedding or projection.
#'
#' @return Code contains the learned shared encoding space. The encoded space refers to the full dimension reduction of both samples and features after matrix multiplication by parameters K and v for y, as well as, L and u for x. The decode is an estimation of the full matrix dataset, where the code is used and matrix multiplied as t(K)(Y_code)t(v), and t(L)(X_code)t(u) to calculate the decoded estimation.
#'
#' @return Recover contains the list of predictions for the test dataset as indicated by a 1 in the binary prediction matrices. The prediction occurs in the shared lower dimensional space where all data sets in the list are projected to using a common latent code.
#'
#' @export
gcode <- function(data_list,
                  config = gcode::extract_config(verbose = F),
                  transfer = gcode::extract_transfer_framework(verbose = F),
                  recover = gcode::extract_recovery_framework(verbose = F),
                  join = gcode::extract_join_framework(verbose=F)
){
  
  runtime.start <- Sys.time()
  
  set.seed(config$seed)
  
  convergence.parameters <- list(count=0,score.vec=c())
  
  recover$predict.list <- lapply(c(1:length(data_list)),function(X){NULL})
  
  initialise.model <- initialise.gcode(data_list = data_list,
                                       config = config,
                                       transfer = transfer,
                                       join = join)
  
  main.parameters <- initialise.model$main.parameters
  main.code <- initialise.model$main.code
  common.template <- initialise.model$common.template
  
  if (config$verbose){
    print(paste("Beginning gcode learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim, "    Latent invariant dimension (config$k_dim): ", config$k_dim, "    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose, sep=""))
  }
  
  
  while (T){
    
    prev_code <- main.code
    
    for (i in 1:length(data_list)){
      
      internal.parameters <- list(alpha=main.parameters$alpha[[join$complete$alpha[i]]],
                                  beta=main.parameters$beta[[join$complete$beta[i]]])
      
      internal.code <- list(encode=main.code$encode[[join$complete$code[i]]],
                            code=main.code$code[[join$complete$code[i]]])
      
      internal.parameters$alpha[,join$labels$alpha] <- common.template$alpha
      internal.parameters$beta[,join$labels$beta] <- common.template$beta
      
      return_update <- update_set(x = as.matrix(data_list[[i]]),
                                  main.parameters = internal.parameters,
                                  main.code = internal.code, 
                                  config = config,
                                  fix = transfer$fix
                                  )
      
      main.parameters$alpha[[join$complete$alpha[i]]] <- return_update$main.parameters$alpha
      main.parameters$beta[[join$complete$beta[i]]] <- return_update$main.parameters$beta
      
      common.template$alpha <- main.parameters$alpha[[join$complete$alpha[i]]][,row.names(data_list[[i]])%in%join$labels$alpha]
      common.template$beta <- main.parameters$beta[[join$complete$beta[i]]][colnames(data_list[[i]])%in%join$labels$beta,]
      
      main.parameters$alpha[[join$complete$alpha[i]]][,row.names(data_list[[i]])%in%join$labels$alpha] <- return_update$main.parameters$alpha[,row.names(data_list[[i]])%in%join$labels$alpha]
      main.parameters$bera[[i]][colnames(data_list[[i]])%in%join$labels$beta,] <- return_update$main.parameters$beta[colnames(data_list[[i]])%in%join$labels$beta,]
      
      main.code$code[[join$complete$code[i]]] <- return_update$main.code$code
      main.code$encode[[join$complete$code[i]]] <- return_update$main.code$encode

    }
    
    matrix.residuals <- Reduce("+",main.code$encode) - Reduce("+",prev_code$encode)
    
    total.mse <- mean(abs(matrix.residuals))
    
    # Check convergence
    convergence.parameters$score.vec <- c(convergence.parameters$score.vec, total.mse)
    MSE <- tail(convergence.parameters$score.vec,1)
    prev.MSE <- tail(convergence.parameters$score.vec,2)[1]
    
    if (convergence.parameters$count>=1){
      print(paste("Iteration: ",convergence.parameters$count," with Tolerance of: ", abs(prev.MSE - MSE),sep=""))
      if ((convergence.parameters$count >= config$max_iter ) | abs(prev.MSE - MSE) < config$tol){
        break
      }
    }
    
    
    convergence.parameters$count = convergence.parameters$count + 1
    
    if (any(do.call('c',lapply(recover$design.list,function(X){!is.null(X)})))){
      
      recover_data <- recover_points(
        data_list,
        main.code = main.code,
        main.parameters = main.parameters,
        config = config,
        recover = recover,
        join = join
      )
      
      recover <- recover_data$recover
      data_list <- recover_data$data_list
      
    }
    
  }
  
  if (config$verbose){
    print("Learning has converged for gcode, beginning (if requested) dimension reduction")
  }
  
  return_list <- list(
    
    main.parameters = main.parameters,
    
    main.code = main.code,
    
    common.template = common.template,
    
    recover =  recover,
    
    meta.parameters = list(
      config = config,
      join = join
    ),
    
    convergence.parameters = convergence.parameters
    
    
  )
  
  if (config$dimension_reduction){
    dimension_reduction <- lapply(c(1:length(data_list)),function(Y){
      
      x <- as.matrix(data_list[[Y]])
      
      feature_x.dim_reduce.encode <- t(main.parameters$alpha[[join$complete$alpha[i]]]%*%x)
      sample_x.dim_reduce.encode <- x%*%main.parameters$beta[[join$complete$beta[i]]]
      
      return(list(
        feature_x.dim_reduce.encode = feature_x.dim_reduce.encode,
        sample_x.dim_reduce.encode = sample_x.dim_reduce.encode
      ))
    })
    
    return_list$dimension_reduction <- dimension_reduction
  }
  
  runtime.end <- Sys.time()
  
  return_list$runtime <- list(runtime.start = runtime.start,
                              runtime.end = runtime.end,
                              runtime.total = runtime.end - runtime.start)
  
  if (config$verbose){
    print(paste("Done! Total runtime of   ", runtime.end - runtime.start ,sep=""))
  }
  
  return(return_list)
  
}

#' @export
update_set <- function(x,
                       main.parameters,
                       main.code, 
                       config,
                       fix){
  
  main.parameters$alpha <- if(fix$alpha){main.parameters$alpha}else{(t((x)%*%t((main.code$code)%*%t(main.parameters$beta))%*%pinv(t((main.code$code)%*%t(main.parameters$beta)))))}
  main.parameters$beta <- if(fix$beta){main.parameters$beta}else{(t(pinv(((t(main.parameters$alpha)%*%(main.code$code))))%*%t(t(main.parameters$alpha)%*%(main.code$code))%*%(x)))}
  
  main.parameters$alpha <- if(fix$alpha){main.parameters$alpha}else{soft_threshold(main.parameters$alpha,config)}
  main.parameters$beta <- if(fix$beta){main.parameters$beta}else{soft_threshold( main.parameters$beta,config)}
  
  main.code$encode <- if(fix$encode){main.code$encode}else{(main.parameters$alpha%*%(x)%*%(main.parameters$beta))}
  main.code$code <- if(fix$code){main.code$code}else{pinv(t(main.parameters$alpha))%*%main.code$encode%*%pinv(main.parameters$beta)}
  
  return(list(main.parameters = main.parameters,
              main.code = main.code
  ))
  
}


pinv <- function(X){
  MASS::ginv(t(X)%*%X)
}


chunk <- function(x,n){
  if (n==1){
    list(x)
  }
  else{
    split(x, cut(seq_along(x), n, labels = FALSE))
  }
}

#' @export
soft_threshold <- function(param,config){
  
  alpha <- config$regularise$a
  lambda <- config$regularise$l
  
  gamma <- lambda*alpha
  
  param <- (((param - gamma)*(param>0)+
               (param + gamma)*(param<0))*(abs(param)>gamma))
  
  return(param/(1+lambda*(1-alpha)))
  
}

