#' Penalized Version of Cox Model Prediction
#' 
#' @import survival
#' @import glmnet
#' @param connectome a array indicating the connectivity between M1 edges and over N subjects. The dimension should be `M1*M1*N`.
#' @param x non-image covariates matrix, of `n (obs)* p (vars)`.
#' @param time the follow-up time for all individuals.
#' @param status the status indicator, normally 0=alive and 1=event.
#' @param cv a character indicating the method of cross-validation. The default method is "leave-one-out" cross-validation.
#' @param k a parameter used to set the number of folds for k-fold cross-validation.
#' @param thresh the value of the threshold for selecting significantly related edges. The default value is .01.
#' @param edge a character indicating the model is fitted with either positive and negative edges respectively or combined edges together. The default is "separate".
#' @param type type of penalty. “lasso” for LASSO, “ridge” for ridge, and “EN” for elastic net. The default is “lasso”.
#' @param lambda a user-specified lambda sequence or the optimal one automatically searched by cv.glmnet.
#' @param alpha the elastic net mixing parameter, ranging from 0 to 1. The default is 0.95.
#' @param seed the value used to set seed for random sampling in the process of cross-validation. The default value is 1220.
#' @return A list contains positive edges, negative edges, predicted survival, actual survival outcomes, and the optimal lambda.
#' @export


cox.regularized.GenCPM <- function(connectome, x=NULL, time, status,
                    cv="leave-one-out", k = dim(connectome)[3], 
                    thresh = .01, edge="separate", type="lasso", 
                    lambda=NULL, alpha=NULL, seed = 1220){
  
  if(length(k)>1){
    stop("invalid cross validation index")
  }
  
  if(type=="lasso"){
    type_alpha <- 1
  }else if(type=="ridge"){
    type_alpha <- 0
  }else if(type=="EN"){
    if(is.null(alpha)){
      type_alpha <- 0.95
    }else {
      type_alpha <- alpha
    }
  }else{
    stop("penalty type should only be lasso, ridge or elastic net")
  }
  
  # parse the dimensions of data
  M1 <- dim(connectome)[1]
  N <- dim(connectome)[3]
  M <- M1*(M1-1)/2
  
  # convert connectome to edge matrix
  
  all_edges <- matrix(NA,nrow=M,ncol=N)
  
  for (i in 1:N){
    all_edges[,i] <- t(connectome[,,i])[lower.tri(connectome[,,i],diag=FALSE)]
  }
  
  if((cv == "leave-one-out" && k == dim(connectome)[3]) | (cv == "k-fold" && k == dim(connectome)[3])){
    
    coef <- rep(NA,M)
    pval <- rep(NA,M)
    
    selected_edges_pos <- vector(mode='list', length=N)
    selected_edges_neg <- vector(mode='list', length=N)

    lambda_total <- rep(0,N)
    lambda_total_pos <- rep(0,N)
    lambda_total_neg <- rep(0,N)
    
    lp_pred_pos <- rep(0,N)
    lp_pred_neg <- rep(0,N)
    lp_pred <- rep(0,N)
    
    status_actual <- rep(0,N)
    time_actual <- rep(0,N)

    for(leftout in 1:N){
      
      # divide into TRAIN set and TEST set
      
      train_mats <- all_edges[,-leftout]
      train_time <- time[-leftout]
      train_status <- status[-leftout]
      
      if(missing(x)){
        train_x <- NULL
      }else if(dim(x)[2] == 1){
        train_x <- x[-leftout]
      }else if(dim(x)[2] > 1){
        train_x <- x[-leftout,]
      }else{
        stop("x should be a DataFrame")
      }
      
      test_mats <- all_edges[,leftout]
      test_time <- time[leftout]
      test_status <- status[leftout]
      
      if(missing(x)){
        test_x <- NULL
      }else if(dim(x)[2] == 1){
        test_x <- x[leftout]
      }else if(dim(x)[2] > 1){
        test_x <- x[leftout,]
      }else{
        stop("x should be a DataFrame")
      }
      
      for(j in 1:M){
        ms_edge <- data.frame(time=train_time, status=train_status, edge=train_mats[j,])
        fit_cox <- coxph(Surv(time,status) ~ ., data=ms_edge)
        coef[j] <- summary(fit_cox)$coefficients[1]
        pval[j] <- summary(fit_cox)$coefficients[5]
      }
      
      pos_edge_index <- which(coef > 0 & pval<=thresh)
      neg_edge_index <- which(coef < 0 & pval<=thresh)
      edge_index <- which(pval<=thresh)
      
      nedge <- length(edge_index)
      nedge_pos <- length(pos_edge_index)
      nedge_neg <- length(neg_edge_index)
      
      if(edge=="separate"){
        
        train_mats_pos <- train_mats[which(coef > 0 & pval<=thresh),]
        train_mats_neg <- train_mats[which(coef < 0 & pval<=thresh),]
        
        test_mats_pos <- test_mats[which(coef > 0 & pval<=thresh)]
        test_mats_neg <- test_mats[which(coef < 0 & pval<=thresh)]
        
        train_covars_pos <- data.matrix(cbind(edge=t(train_mats_pos), x=train_x))
        train_covars_neg <- data.matrix(cbind(edge=t(train_mats_neg), x=train_x))
        
        test_covars_pos <- data.matrix(cbind(edge=t(test_mats_pos), x=test_x))      
        test_covars_neg <- data.matrix(cbind(edge=t(test_mats_neg), x=test_x))      
        
        if(missing(lambda) | length(lambda) < 1){
          
          fit_lambda_pos <- cv.glmnet(train_covars_pos, Surv(train_time, train_status), 
                                  family="cox", alpha=type_alpha, 
                                  penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos))) # feature selection is only for edges, not for non-image variables
          
          opt_lambda_pos <- fit_lambda_pos$lambda.1se
          
          fit_lambda_neg <- cv.glmnet(train_covars_neg, Surv(train_time, train_status), 
                                      family="cox", alpha=type_alpha, 
                                      penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
          
          opt_lambda_neg <- fit_lambda_neg$lambda.1se
          
        }else if(length(lambda)>1){
          
          fit_lambda_pos <- cv.glmnet(train_covars_pos, Surv(train_time, train_status), 
                                  family="cox", lambda = lambda, 
                                  alpha=type_alpha,
                                  penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
          
          opt_lambda_pos <- fit_lambda_pos$lambda.1se
          
          fit_lambda_neg <- cv.glmnet(train_covars_neg, Surv(train_time, train_status), 
                                      family="cox", lambda = lambda, 
                                      alpha=type_alpha,
                                      penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
          
          opt_lambda_neg <- fit_lambda_neg$lambda.1se
          
        }else{
          
          opt_lambda_pos <- lambda
          opt_lambda_neg <- lambda
          
        }
        
        lambda_total_pos[leftout] <- opt_lambda_pos
        lambda_total_neg[leftout] <- opt_lambda_neg
        
        fit_pos <- glmnet(train_covars_pos, Surv(train_time, train_status), 
                      family="cox", lambda = opt_lambda_pos, alpha=type_alpha,
                      penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
        
        fit_neg <- glmnet(train_covars_neg, Surv(train_time, train_status), 
                          family="cox", lambda = opt_lambda_neg, alpha=type_alpha,
                          penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
        
        # edge selection by glmnet
        
        coef_edge_pos <- data.matrix(coef(fit_pos))[c(1:nedge_pos),]
        coef_edge_neg <- data.matrix(coef(fit_neg))[c(1:nedge_neg),]
        
        selected_edges_pos[[leftout]] <- pos_edge_index[which(coef_edge_pos > 0)]
        selected_edges_neg[[leftout]] <- neg_edge_index[which(coef_edge_neg < 0)]
        
        # predict TEST sub survival
        
        lp_pred_pos[leftout] <- predict(fit_pos, newx=test_covars_pos, s = opt_lambda_pos)
        lp_pred_neg[leftout] <- predict(fit_neg, newx=test_covars_neg, s = opt_lambda_neg)
        
        status_actual[leftout] <- test_status
        time_actual[leftout] <- test_time
        
        
      }else if(edge=="combined"){
        
        train_mats <- train_mats[which(pval<=thresh),]
        test_mats <- test_mats[which(pval<=thresh)]
        
        train_covars <- data.matrix(cbind(edge=t(train_mats), x=train_x))
        test_covars <- data.matrix(cbind(edge=t(test_mats), x=test_x))      
        
        if(missing(lambda) | length(lambda) < 1){
          
          fit_lambda <- cv.glmnet(train_covars, Surv(train_time, train_status), 
                                  family="cox", alpha=type_alpha,
                                  penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
          opt_lambda <- fit_lambda$lambda.1se
          
        }else if(length(lambda)>1){
          
          fit_lambda <- cv.glmnet(train_covars, Surv(train_time, train_status), 
                                  family="cox", lambda = lambda, 
                                  alpha=type_alpha,
                                  penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
          opt_lambda <- fit_lambda$lambda.1se
          
        }else{
          
          opt_lambda <- lambda
          
        }
        
        lambda_total[leftout] <- opt_lambda
        
        fit <- glmnet(train_covars, Surv(train_time, train_status), 
                      family="cox", lambda = opt_lambda, alpha=type_alpha,
                      penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
        
        # edge selection by glmnet
        
        coef_edge <- data.matrix(coef(fit))[c(1:nedge),]
        
        selected_edges_pos[[leftout]] <- edge_index[which(coef_edge > 0)]
        selected_edges_neg[[leftout]] <- edge_index[which(coef_edge < 0)]
        
        # predict TEST sub survival
        
        lp_pred[leftout] <- predict(fit, newx=test_covars, s = opt_lambda)
        
        status_actual[leftout] <- test_status
        time_actual[leftout] <- test_time

      }
      
    }
    
    if(edge=="separate"){
      
      return(list(positive_edges=selected_edges_pos, 
                  negative_edges=selected_edges_neg,
                  positive_predicted_linear_predictor=lp_pred_pos,
                  negative_predicted_linear_predictor=lp_pred_neg,
                  actual_status=status_actual, actual_time=time_actual,
                  positive_lambda_total=lambda_total_pos, 
                  negative_lambda_total=lambda_total_neg))
      
    }else if(edge=="combined"){
      
      return(list(positive_edges=selected_edges_pos, 
                  negative_edges=selected_edges_neg,
                  predicted_linear_predictor=lp_pred, 
                  actual_status=status_actual, actual_time=time_actual,
                  lambda_total=lambda_total))
    }
    
  }
  
  if(cv == "k-fold" | ((cv == "leave-one-out") && (k != dim(connectome)[3]))){
    
    if(k-floor(k) != 0){
      stop("k should be an integer.")
    }
    else{
      
      set.seed(seed)
      
      # stratified sampling for TRAIN and TEST subs
      
      r0 <- length(which(status==0))/N
      r1 <- length(which(status==1))/N
      
      n0 <- floor(N/k * r0)
      n1 <- floor(N/k * r1)
      
      samplesize <- n0+n1
      
      index_list <- vector(mode="list",length=2)
      final_index <- matrix(NA,nrow=k,ncol=samplesize)
      
      index_list[[1]] <- which(status==0)
      index_list[[2]] <- which(status==1)
      
      for (fold in 1:k){
        
        randinds0 <- sample(index_list[[1]], size=n0, replace=FALSE)
        randinds1 <- sample(index_list[[2]], size=n1, replace=FALSE)
        
        final_index[fold,] <- sample(c(randinds0, randinds1),size=samplesize)
        
        index_list[[1]] <- setdiff(index_list[[1]],randinds0)
        index_list[[2]] <- setdiff(index_list[[2]],randinds1)
        
      }
      
      coef <- rep(NA,M)
      pval <- rep(NA,M)
      
      selected_edges_pos <- vector(mode='list', length=k)
      selected_edges_neg <- vector(mode='list', length=k)
      
      lambda_total <- rep(0,k)
      lambda_total_pos <- rep(0,k)
      lambda_total_neg <- rep(0,k)
      
      lp_pred <- matrix(0,k,samplesize)
      lp_pred_pos <- matrix(0,k,samplesize)
      lp_pred_neg <- matrix(0,k,samplesize)
      
      status_actual <- matrix(0,k,samplesize)
      time_actual <- matrix(0,k,samplesize)

      for (fold in 1:k){
        
        test_index <- final_index[fold,]
        train_index <- as.vector(final_index[-fold,])
        
        # divide into TRAIN set and TEST set
        
        train_mats <- all_edges[,train_index]
        train_time <- time[train_index]
        train_status <- status[train_index]
        
        if(missing(x)){
          train_x <- NULL
        }else if(dim(x)[2] == 1){
          train_x <- x[train_index]
        }else if(dim(x)[2] > 1){
          train_x <- x[train_index,]
        }else{
          stop("x should be a DataFrame")
        }
        
        test_mats <- all_edges[,test_index]
        test_time <- time[test_index]
        test_status <- status[test_index]
        
        if(missing(x)){
          test_x <- NULL
        }else if(dim(x)[2] == 1){
          test_x <- x[test_index]
        }else if(dim(x)[2] > 1){
          test_x <- x[test_index,]
        }else{
          stop("x should be a DataFrame")
        }
        
        # marginal screening
        
        for(j in 1:M){
          ms_edge <- data.frame(time=train_time, status=train_status, edge=train_mats[j,])
          fit_cox <- coxph(Surv(time,status) ~ ., data=ms_edge)
          coef[j] <- summary(fit_cox)$coefficients[1]
          pval[j] <- summary(fit_cox)$coefficients[5]
        }
        
        edge_index <- which(pval<=thresh)
        pos_edge_index <- which(coef > 0 & pval<=thresh)
        neg_edge_index <- which(coef < 0 & pval<=thresh)
        
        nedge <- length(edge_index)
        nedge_pos <- length(pos_edge_index)
        nedge_neg <- length(neg_edge_index)
        
        if(edge=="separate"){
          
          train_mats_pos <- train_mats[which(coef > 0 & pval<=thresh),]
          train_mats_neg <- train_mats[which(coef < 0 & pval<=thresh),]
          
          test_mats_pos <- test_mats[which(coef > 0 & pval<=thresh),]
          test_mats_neg <- test_mats[which(coef < 0 & pval<=thresh),]
          
          train_covars_pos <- data.matrix(cbind(edge=t(train_mats_pos), x=train_x))
          train_covars_neg <- data.matrix(cbind(edge=t(train_mats_neg), x=train_x))
          
          test_covars_pos <- data.matrix(cbind(edge=t(test_mats_pos), x=test_x))
          test_covars_neg <- data.matrix(cbind(edge=t(test_mats_neg), x=test_x))
          
          if(missing(lambda) | length(lambda) < 1){
            
            fit_lambda_pos <- cv.glmnet(train_covars_pos, Surv(train_time, train_status), 
                                    family="cox", alpha=type_alpha, 
                                    penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
            opt_lambda_pos <- fit_lambda_pos$lambda.1se
            
            fit_lambda_neg <- cv.glmnet(train_covars_neg, Surv(train_time, train_status), 
                                        family="cox", alpha=type_alpha, 
                                        penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
            opt_lambda_neg <- fit_lambda_neg$lambda.1se
            
          }else if(length(lambda)>1){
            
            fit_lambda_pos <- cv.glmnet(train_covars_pos, Surv(train_time, train_status), 
                                    family="cox", lambda = lambda, 
                                    alpha=type_alpha, 
                                    penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
            opt_lambda_pos <- fit_lambda_pos$lambda.1se
            
            fit_lambda_neg <- cv.glmnet(train_covars_neg, Surv(train_time, train_status), 
                                        family="cox", lambda = lambda, 
                                        alpha=type_alpha,
                                        penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
            opt_lambda_neg <- fit_lambda_neg$lambda.1se
            
          }else{
            
            opt_lambda_pos <- lambda
            opt_lambda_neg <- lambda
            
          }
          
          lambda_total_pos[fold] <- opt_lambda_pos
          lambda_total_neg[fold] <- opt_lambda_neg
          
          fit_pos <- glmnet(train_covars_pos, Surv(train_time, train_status), 
                        family="cox", lambda = opt_lambda_pos, alpha=type_alpha,
                        penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
          
          fit_neg <- glmnet(train_covars_neg, Surv(train_time, train_status), 
                            family="cox", lambda = opt_lambda_neg, 
                            alpha=type_alpha,
                            penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
          
          # edge selection by glmnet
          
          coef_edge_pos <- data.matrix(coef(fit_pos))[c(1:nedge_pos),]
          coef_edge_neg <- data.matrix(coef(fit_neg))[c(1:nedge_neg),]
          
          selected_edges_pos[[fold]] <- pos_edge_index[which(coef_edge_pos > 0)]
          selected_edges_neg[[fold]] <- neg_edge_index[which(coef_edge_neg < 0)]
          
          # predict TEST sub survival
          
          lp_pred_pos[fold,] <- predict(fit_pos, newx=test_covars_pos, s=opt_lambda_pos)
          lp_pred_neg[fold,] <- predict(fit_neg, newx=test_covars_neg, s=opt_lambda_pos)
          
          status_actual[fold,] <- test_status
          time_actual[fold,] <- test_time
          
          
        }else if(edge=="combined"){
          
          train_mats <- train_mats[which(pval<=thresh),]
          test_mats <- test_mats[which(pval<=thresh),]
          
          train_covars <- data.matrix(cbind(edge=t(train_mats), x=train_x))
          test_covars <- data.matrix(cbind(edge=t(test_mats), x=test_x))
          
          if(missing(lambda) | length(lambda) < 1){
            
            fit_lambda <- cv.glmnet(train_covars, Surv(train_time, train_status), 
                                    family="cox", alpha=type_alpha, 
                                    penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
            opt_lambda <- fit_lambda$lambda.1se
            
          }else if(length(lambda)>1){
            
            fit_lambda <- cv.glmnet(train_covars, Surv(train_time, train_status), 
                                    family="cox", lambda = lambda, 
                                    alpha=type_alpha, 
                                    penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
            opt_lambda <- fit_lambda$lambda.1se
            
          }else{
            
            opt_lambda <- lambda
            
          }
          
          lambda_total[fold] <- opt_lambda
          
          fit <- glmnet(train_covars, Surv(train_time, train_status), 
                        family="cox", lambda = opt_lambda, alpha=type_alpha,
                        penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
          
          # edge selection by glmnet
          
          coef_edge <- data.matrix(coef(fit))[c(1:nedge),]

          selected_edges_pos[[fold]] <- edge_index[which(coef_edge > 0)]
          selected_edges_neg[[fold]] <- edge_index[which(coef_edge < 0)]

          # predict TEST sub survival
          
          lp_pred[fold,] <- predict(fit, newx=test_covars, s=opt_lambda)
          status_actual[fold,] <- test_status
          time_actual[fold,] <- test_time
          
        }
        
      }
      
      if(edge=="separate"){
        
        return(list(positive_edges=selected_edges_pos, 
                    negative_edges=selected_edges_neg,
                    positive_predicted_linear_predictor=lp_pred_pos,
                    negative_predicted_linear_predictor=lp_pred_neg,
                    actual_status=status_actual, actual_time=time_actual,
                    positive_lambda_total=lambda_total_pos, 
                    negative_lambda_total=lambda_total_neg))
        
      }else if(edge=="combined"){
        
        return(list(positive_edges=selected_edges_pos, 
                    negative_edges=selected_edges_neg,
                    predicted_linear_predictor=lp_pred, 
                    actual_status=status_actual, actual_time=time_actual,
                    lambda_total=lambda_total))
      }
      
    }
    
    
  }
  
}

