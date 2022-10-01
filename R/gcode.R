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
                  join = gcode::extract_join_framework(verbose=F),
                  references = gcode::extract_references_framework(F)
){

  runtime.start <- Sys.time()

  set.seed(config$seed)

  convergence.parameters <- list(count=0,score.vec=c())

  initialise.model <- initialise.gcode(data_list = data_list,
                                       config = config,
                                       transfer = transfer,
                                       join = join)

  main.parameters <- initialise.model$main.parameters
  main.code <- initialise.model$main.code

  
  if (config$verbose){
    print(paste("Beginning gcode learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim, "    Latent invariant dimension (config$k_dim): ", config$k_dim, "    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose, sep=""))
  }

  
  while (T){

    prev.code <- main.code
    
    for (i in 1:length(join$complete$data_list)){
      
      
      internal.parameters <- list(alpha_sample=main.parameters$alpha_sample[[join$complete$alpha_sample[i]]],
                                  beta_sample=main.parameters$beta_sample[[join$complete$beta_sample[i]]]
      )
      
      internal.code <- list(encode=main.code$encode[[join$complete$code[i]]],
                            code=main.code$code[[join$complete$code[i]]])
      
      reference.parameters <- list(alpha_sample=main.parameters$alpha_sample[[join$complete$alpha_sample[references$data_list]]],
                                   beta_sample=main.parameters$beta_sample[[join$complete$beta_sample[references$data_list]]]
      )
      
      reference.code <- list(encode=main.code$encode[[join$complete$code[references$data_list]]],
                             code=main.code$code[[join$complete$code[references$data_list]]])
      
      return_update <- update_regression(y = data_list[[references$data_list]],
                                         x = as.matrix(data_list[[join$complete$data_list[i]]]),
                                         main.parameters = internal.parameters,
                                         main.code = internal.code,
                                         reference.parameters = reference.parameters,
                                         reference.code = reference.code,
                                         config = config,
                                         fix = transfer$fix,
                                         seed_iter = convergence.parameters$count
      )
      
      main.parameters$alpha_sample[[join$complete$alpha_sample[i]]] <- return_update$main.parameters$alpha_sample
      main.parameters$beta_sample[[join$complete$beta_sample[i]]] <- if(join$covariance[i]){t(return_update$main.parameters$alpha_sample)}else{return_update$main.parameters$beta_sample}
      
      
    }
    

    for (i in 1:length(join$complete$data_list)){

      internal.parameters <- list(alpha_sample=main.parameters$alpha_sample[[join$complete$alpha_sample[i]]],
                                  beta_sample=main.parameters$beta_sample[[join$complete$beta_sample[i]]]
                                  )
      
      internal.code <- list(encode=main.code$encode[[join$complete$code[i]]],
                            code=main.code$code[[join$complete$code[i]]])
      
      return_update <- update_set(x = as.matrix(data_list[[join$complete$data_list[i]]]),
                                  main.parameters = internal.parameters,
                                  main.code = internal.code,
                                  config = config,
                                  fix = transfer$fix,
                                  seed_iter = convergence.parameters$count
      )
      
      main.parameters$alpha_sample[[join$complete$alpha_sample[i]]] <- return_update$main.parameters$alpha_sample
      main.parameters$beta_sample[[join$complete$beta_sample[i]]] <- if(join$covariance[i]){t(return_update$main.parameters$alpha_sample)}else{return_update$main.parameters$beta_sample}
      
      main.code$code[[join$complete$code[i]]] <- return_update$main.code$code
      main.code$encode[[join$complete$code[i]]] <- return_update$main.code$encode
      
    }
    
    
    
    

    total.mse <- matrix.residuals <- Reduce("+",lapply(unique(join$complete$code),function(X){

      mean(abs(main.code$encode[[X]] - prev.code$encode[[X]]))

    }))/length(join$complete$code)

    # Check convergence
    convergence.parameters$score.vec <- c(convergence.parameters$score.vec, total.mse)
    MSE <- mean(tail(convergence.parameters$score.vec,1))
    prev.MSE <- mean(tail(convergence.parameters$score.vec,2)[1])
    
    if (convergence.parameters$count>=1){
      if (config$verbose){
        print(paste("Iteration: ",convergence.parameters$count," with Tolerance of: ", abs(prev.MSE - MSE),sep=""))
      }
      
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
        join = join,
        references = references,
        iter = convergence.parameters$count
      )
      
      if (config$verbose){
        print(c("Recovery operation round:    ", convergence.parameters$count))
      }
      
      
      recover$predict.list <- data_list <- recover_data$data_list
      
    }
    

  }

  if (config$verbose){
    print("Learning has converged for gcode, beginning (if requested) dimension reduction")
  }
  
  
  recover$predict.list <- lapply(c(1:length(data_list)),function(X){
    Matrix::Matrix(recover$design.list[[X]]*data_list[[join$complete$data_list[X]]],sparse=T)
  })

  return_list <- list(

    main.parameters = main.parameters,

    main.code = main.code,

    recover =  recover,

    meta.parameters = list(
      config = config,
      join = join
    ),

    convergence.parameters = convergence.parameters


  )

  if (config$dimension_reduction){
    dimension_reduction <- lapply(c(1:length(join$complete$data_list)),function(Y){

      x <- as.matrix(data_list[[join$complete$data_list[Y]]])

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
                       fix,
                       seed_iter){

  set.seed(seed_iter)
  chunk_list <- list(alpha=chunk(sample(c(1:config$i_dim)),round(config$i_dim/config$batch_size,0)),
                     beta=chunk(sample(c(1:config$j_dim)),round(config$j_dim/config$batch_size,0)))
  
  internal_parameters <- main.parameters
  internal_code <- main.code
  
  param_combin <- cbind(rep(1:round(config$i_dim/config$batch_size,0),each=round(config$j_dim/config$batch_size,0)),rep(1:round(config$j_dim/config$batch_size,0),times=round(config$i_dim/config$batch_size,0)))

  return_update <- parallel::mclapply(c(1:dim(param_combin)[1]),function(i){
    
    beta_sample <- if(fix$beta_sample){main.parameters$beta_sample[,chunk_list$beta[[param_combin[i,2]]]]}else{
      (1-config$learn_rate)*
        internal_parameters$beta_sample[,chunk_list$beta[[param_combin[i,2]]],drop=F]+
        config$learn_rate*
        array(t(MASS::ginv(t((internal_code$code[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]],drop=F]))%*%
                                 ((internal_code$code[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]],drop=F])))%*%
                      t(internal_code$code[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]],drop=F])%*%
                      MASS::ginv(internal_parameters$alpha_sample[chunk_list$alpha[[param_combin[i,1]]],,drop=F]%*%
                                   t(internal_parameters$alpha_sample[chunk_list$alpha[[param_combin[i,1]]],,drop=F]))%*%
                      internal_parameters$alpha_sample[chunk_list$alpha[[param_combin[i,1]]],,drop=F]%*%
                      (x)),dim=c(dim(x)[2],length(chunk_list$beta[[param_combin[i,2]]])))}
    alpha_sample <- if(fix$alpha_sample){main.parameters$alpha_sample[chunk_list$alpha[[param_combin[i,1]]],]}else{
      (1-config$learn_rate)*
        internal_parameters$alpha_sample[chunk_list$alpha[[param_combin[i,1]]],]+
        config$learn_rate*
        array(t((x)%*%
                      internal_parameters$beta_sample[,chunk_list$beta[[param_combin[i,2]]],drop=F]%*%
                      MASS::ginv(t(internal_parameters$beta_sample[,chunk_list$beta[[param_combin[i,2]]],drop=F])%*%
                                   internal_parameters$beta_sample[,chunk_list$beta[[param_combin[i,2]]],drop=F])%*%
                      t(internal_code$code[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]],drop=F])%*%
                      MASS::ginv((internal_code$code[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]]])%*%
                                   t(internal_code$code[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]]]))),dim = c(length(chunk_list$alpha[[param_combin[i,1]]]),dim(x)[1]))}
    
    
    encode <- if(fix$encode){main.code$encode}else{(internal_parameters$alpha_sample[chunk_list$alpha[[param_combin[i,1]]],]%*%(x)%*%(internal_parameters$beta_sample[,chunk_list$beta[[param_combin[i,2]]]]))}

    code <- if(fix$code){main.code$code[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]]]}else{
      (1-config$learn_rate)*
        internal_code$code[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]]]+
        config$learn_rate*
        MASS::ginv(internal_parameters$alpha_sample[chunk_list$alpha[[param_combin[i,1]]],,drop=F]%*%
                     t(internal_parameters$alpha_sample[chunk_list$alpha[[param_combin[i,1]]],,drop=F]))%*%
        internal_code$encode[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]]]%*%
        MASS::ginv(t(internal_parameters$beta_sample[,chunk_list$beta[[param_combin[i,2]]],drop=F])%*%
                     internal_parameters$beta_sample[,chunk_list$beta[[param_combin[i,2]]],drop=F])}

    return(list(main.parameters = list(alpha_sample=alpha_sample,
                                       beta_sample=beta_sample),
                main.code = list(
                  encode=encode,
                  code=code
                )
    ))
    
  },mc.cores = config$n.cores)
  
  for (i in c(1:dim(param_combin)[1])){
    
    main.parameters$alpha_sample[chunk_list$alpha[[param_combin[i,1]]],] <- if(fix$alpha_sample){main.parameters$alpha_sample[chunk_list$alpha[[param_combin[i,1]]],]}else{soft_threshold(return_update[[i]]$main.parameters$alpha_sample,config = config)}
    main.parameters$beta_sample[,chunk_list$beta[[param_combin[i,2]]]] <- if(fix$beta_sample){main.parameters$beta_sample[,chunk_list$beta[[param_combin[i,2]]]]}else{soft_threshold(return_update[[i]]$main.parameters$beta_sample,config = config)}

    main.code$encode[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]]] <- if(fix$code){main.code$encode[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]]]}else{return_update[[i]]$main.code$encode}
    main.code$code[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]]] <- if(fix$code){main.code$code[chunk_list$alpha[[param_combin[i,1]]],chunk_list$beta[[param_combin[i,2]]]]}else{return_update[[i]]$main.code$code}
    
  }

  return(list(main.parameters = main.parameters,
              main.code = main.code
  ))

}






update_regression <- function(y,
                              x,
                              main.parameters,
                              main.code,
                              reference.parameters,
                              reference.code,
                              config,
                              fix,
                              seed_iter){
  
  for (i in 1:3){
  reference.encode <- (main.parameters$alpha_sample)%*%y%*%(reference.parameters$beta_sample)
  
  main.encode <- (main.parameters$alpha_sample)%*%x%*%(main.parameters$beta_sample)
  D <- MASS::ginv(t(main.encode)%*%(main.encode))%*%t(main.encode)%*%reference.encode
  main.parameters$beta_sample <- main.parameters$beta_sample%*%D
  
  main.encode <- (main.parameters$alpha_sample)%*%x%*%(main.parameters$beta_sample)
  H <- reference.encode%*%t(main.encode)%*%MASS::ginv(main.encode%*%t(main.encode))
  main.parameters$alpha_sample <- H%*%main.parameters$alpha_sample
  }
  return(list(main.parameters = main.parameters))

}

#' @export
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
  
  return(param / (1 + lambda*(1 - alpha)))
  
}