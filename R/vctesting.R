######################################################################
## Editor: Shiquan Sun and Jiaqiang Zhu
## Date: 2018-10-28 22:03:36
## Modified: 2019-4-26 11:24:27
## Affiliation: University of Michigan
######################################################################


#' Testing multiple kernel matrices
#' 
#' @param object SPRINT object
#' @param check_positive Check the kernel matrix is positive or not
#' @param verbose Output fitting information
#' @export
sprint.test <- function(object, check_positive = TRUE, verbose = TRUE) {
  
    # Euclid distance, and compute the range of parameter l
    ED <- as.matrix(dist(object@location[ ,1:2]))
    lrang <- ComputeGaussianPL(ED, compute_distance=FALSE)[3:7]
    
    # check one kernel at a time
    res_kernel <- list()
    res_pval <- NULL
    for(ikernel in c(1:5) ){
      # Gaussian kernel
      kernel_mat <- exp(-ED^2/(2*lrang[ikernel]^2))
      object <- sprint.test_each(object, kernel_mat=kernel_mat, check_positive=check_positive, verbose=verbose)
      res_pval <- cbind(res_pval, object@res_stest$sw)
      res_kernel <- c(res_kernel, object@res_stest)
      rm(kernel_mat)

      # Periodic kernel
      kernel_mat <- cos(2*pi*ED/lrang[ikernel])
      object <- sprint.test_each(object, kernel_mat=kernel_mat, check_positive=check_positive, verbose=verbose)
      res_pval <- cbind(res_pval, object@res_stest$sw)
      res_kernel <- c(res_kernel, object@res_stest)
      rm(kernel_mat)
    }# end for ikernel
    
    colnames(res_pval) <- paste0(c("GSP","COS"), rep(1:5,each=2))
    rownames(res_pval) <- rownames(object@counts)
    
    object@res_stest <- res_kernel
    ## summarize ten pvalues into one
    combined_pvalue <- CombinePValues(res_pval)
    object@res_mtest <- data.frame(res_pval, combined_pvalue = combined_pvalue,  adjusted_pvalue = p.adjust(combined_pvalue, method="BY") )
   
    return(object)
}# end function score test



#' Testing one kernel matrix to identify spatial pattern
#' 
#' @param object SPRINT object
#' @param kernel_mat The kernel matrix
#' @param check_positive Check the kernel matrix is positive or not
#' @param verbose Output fitting information
#' @export
#' 
sprint.test_each <- function(object, kernel_mat, check_positive = FALSE, verbose = TRUE) {
    
  # The number of core used in testing step
  num_core <- object@num_core
  
  if(check_positive){# kernel matrix should be positive definition matrix
    # need to check positive definition before tesing
    eig <- eigen(kernel_mat)
    eigval <- eig$value
    eigvector <- eig$vectors
    if(any(eigval<1e-8)){ 
      #warning("SPRINT.TEST::the kernel matrix is singular, it has been modified!")
      eigval[eigval<1e-8] <- 1e-8
      kernel_mat <- eigvector%*%diag(eigval)%*%t(eigvector)
      # kernel_mat <- as.matrix(nearPD(kernel_mat,corr=T,maxit=500)$mat)
      # kernel_mat <- as.matrix(cov2cor(nearPD(kernel_mat,maxit=500)$mat))
      # kernel_mat <- as.matrix(nearPD(kernel_mat,maxit=500)$mat)
      # kernel_mat <- cov2cor(kernel_mat)
    }# end fi
    rm(eig)
    rm(eigval)
    rm(eigvector)
  }# end fi check
  
  # number of genes to testing
  num_gene_test <- length(object@res_vc)
  #========================================
  # using parallel to test variance components
  registerDoParallel(cores = num_core)
  res_test <-foreach(ig = 1:num_gene_test, .combine=rbind)%dopar%{
    if(verbose) {cat(paste("NO. Gene = ",ig,"\n"))}
    model1 <- object@res_vc[[ig]]
    
    davies_sw <-NA
    converged <- FALSE
    if((class(model1) != "try-error") & (!any(is.na(model1$Y)))){
      converged <- model1$converged
      # extract the analyzed cells
      cov_mat <- kernel_mat[model1$analy_indx, model1$analy_indx]
      if(ncol(model1$X)>1){
        rest <- ComputeTestQuantRcpp_cov(model1$Y, model1$Py, model1$X, cov_mat, model1$D^2, model1$theta)
      }else{
        rest <- ComputeTestQuantRcpp_nocov(model1$Y, model1$Py, cov_mat, model1$D^2, model1$theta)
      }# end fi
      
      #scaledchisq.pvalue 	<- pchisq(rest$S0/rest$kk, rest$df, lower.tail = F)
      ##--------------------
      ## davies_sw
      newInfoM 	<- rest$infoMp1
      # calculate the scale parameter and degree of freedom
      ee_sw <- rest$ee
      kk_sw <- newInfoM/(2*ee_sw)
      df_sw <- 2*ee_sw^2/(newInfoM)
      davies_sw	<- pchisq(rest$S0/kk_sw, df_sw, lower.tail = F)
      if(verbose) {cat(paste("SPRINT.SCORE::SW pvalue 1 = ", davies_sw,"\n"))}
    } # end fi
    # to return
    #return( data.frame(geneid = names(object@res_vc)[ig], p_score = scaledchisq.pvalue, p_davies = p_davies, sw=davies_sw, converged=converged) )
    return( data.frame(geneid = names(object@res_vc)[ig], sw=davies_sw, converged=converged) )
  }# end parallel foreach
  
  #######
  object@res_stest <- res_test
  rm(res_test)
  return(object)
}# end function score test for each kernel

#########################################
#             CODE END                  #
#########################################




