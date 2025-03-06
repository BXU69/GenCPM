#' Fit a Cox Model for Survival Outcome using Connectome-based Predictive Modeling
#' 
#' @import survival
#' @param connectome a array indicating the connectivity between M1 edges and over N subjects. The dimension must be `M1*M1*N`.
#' @param x non-image covariates matrix of `n (obs)* p (vars)`.
#' @param time the follow-up time for all individuals.
#' @param status the status indicator, normally 0=alive and 1=event.
#' @param cv a character indicating the method of cross-validation. The default method is "leave-one-out" cross-validation.
#' @param k a parameter used to set the number of folds for k-fold cross-validation.
#' @param thresh the value of the threshold for selecting significantly related edges. The default value is .01.
#' @param edge a character indicating the model is fitted with either positive and negative edges respectively or combined edges together. The default is "separate".
#' @param seed the value used to set seed for random sampling in the process of cross-validation. The default value is 1220.
#' @return A list contains positive edges, negative edges, predicted survival, and actual survival outcome.
#' @export


cox.GenCPM <- function(connectome, x=NULL, time, status,
                    cv="leave-one-out", k = dim(connectome)[3], 
                    thresh = .01, edge="separate", seed = 1220){
  
  if(length(k)>1){
    stop("invalid cross validation index")
  }
  
  # parse the dimensions of data
  M1 <- dim(connectome)[1]
  N <- dim(connectome)[3]
  M <- M1*(M1-1)/2
  
  all_edges <- matrix(NA,nrow=N,ncol=M)
  
  for (i in 1:N){
    all_edges[i,] <- as.vector(connectome[,,i][upper.tri(connectome[,,i])])
  }
    
  if((cv == "leave-one-out" && k == dim(connectome)[3]) | (cv == "k-fold" && k == dim(connectome)[3])){
    
    coef <- rep(NA,M)
    pval <- rep(NA,M)
    
    pos_edge_index <- vector(mode = "list", length=N)
    neg_edge_index <- vector(mode = "list", length=N)
    edge_index <- vector(mode = "list", length=N)
    
    lp_pred_pos <- rep(0,N)
    lp_pred_neg <- rep(0,N)
    lp_pred <- rep(0,N)
    status_actual <- rep(0,N)
    time_actual <- rep(0,N)
    cindex <- rep(0,N)

    for(leftout in 1:N){
        
      # divide into TRAIN set and TEST set
        
      train_mats <- all_edges[-leftout,]
      train_time <- time[-leftout]
      train_status <- status[-leftout]
        
      test_mats <- all_edges[leftout,]
      test_time <- time[leftout]
      test_status <- status[leftout]
      
      # marginal screening on TRAIN subs

      for(j in 1:M){
        ms_edge <- data.frame(time=train_time, status=train_status, edge=train_mats[,j])
        fit_cox <- coxph(Surv(time,status) ~ ., data=ms_edge)
        coef[j] <- summary(fit_cox)$coefficients[1]
        pval[j] <- summary(fit_cox)$coefficients[5]
      }
      
      pos_edge_index[[leftout]] <- which(coef > 0 & pval<=thresh)
      neg_edge_index[[leftout]] <- which(coef < 0 & pval<=thresh)

      if (edge=="separate"){
        
        train_mats_pos <- rowSums(train_mats[,which(coef > 0 & pval<=thresh)])
        train_mats_neg <- rowSums(train_mats[,which(coef < 0 & pval<=thresh)])
        
        test_mats_pos <- sum(test_mats[which(coef > 0 & pval<=thresh)])
        test_mats_neg <- sum(test_mats[which(coef < 0 & pval<=thresh)])
        
        if(missing(x)){
          train_x <- NULL
          train_df_pos <- data.frame(time=train_time, status=train_status, edge=train_mats_pos)
          train_df_neg <- data.frame(time=train_time, status=train_status, edge=train_mats_neg)
        }else if(dim(x)[2] == 1){
          train_x <- x[-leftout]
          train_df_pos <- data.frame(time=train_time, status=train_status, edge=train_mats_pos, x=train_x)
          train_df_neg <- data.frame(time=train_time, status=train_status, edge=train_mats_neg, x=train_x)
        }else if(dim(x)[2] > 1){
          train_x <- x[-leftout,]
          train_df_pos <- data.frame(time=train_time, status=train_status, edge=train_mats_pos, x=train_x)
          train_df_neg <- data.frame(time=train_time, status=train_status, edge=train_mats_neg, x=train_x)
        }else{
          stop("x should be a DataFrame")
        }
        
        if(missing(x)){
          test_x <- NULL
          test_df_pos <- data.frame(time=test_time, status=test_status, edge=t(test_mats_pos))
          test_df_neg <- data.frame(time=test_time, status=test_status, edge=t(test_mats_neg))
        }else if(dim(x)[2] == 1){
          test_x <- x[leftout]
          test_df_pos <- data.frame(time=test_time, status=test_status, edge=t(test_mats_pos), x=test_x)
          test_df_neg <- data.frame(time=test_time, status=test_status, edge=t(test_mats_neg), x=test_x)
        }else if(dim(x)[2] > 1){
          test_x <- x[leftout,]
          test_df_pos <- data.frame(time=test_time, status=test_status, edge=t(test_mats_pos), x=test_x)
          test_df_neg <- data.frame(time=test_time, status=test_status, edge=t(test_mats_neg), x=test_x)
        }else{
          stop("x should be a DataFrame")
        }
        
        fit_pos <- coxph(Surv(time,status) ~ ., data=train_df_pos)
        fit_neg <- coxph(Surv(time,status) ~ ., data=train_df_neg)
        
        # predict TEST sub survival
        
        lp_pred_pos[leftout] <- predict(fit_pos, newdata=test_df_pos)
        lp_pred_neg[leftout] <- predict(fit_neg, newdata=test_df_neg)
        
        status_actual[leftout] <- test_status
        time_actual[leftout] <- test_time
      
      }else if(edge=="combined"){
      
        train_mats <- rowSums(train_mats[,which(pval<=thresh)])
        test_mats <- sum(test_mats[which(pval<=thresh)])
        
        if(missing(x)){
          train_x <- NULL
          train_df <- data.frame(time=train_time, status=train_status, edge=train_mats)
        }else if(dim(x)[2] == 1){
          train_x <- x[-leftout]
          train_df <- data.frame(time=train_time, status=train_status, edge=train_mats, x=train_x)
        }else if(dim(x)[2] > 1){
          train_x <- x[-leftout,]
          train_df <- data.frame(time=train_time, status=train_status, edge=train_mats, x=train_x)
        }else{
          stop("x should be a DataFrame")
        }
        
        if(missing(x)){
          test_x <- NULL
          test_df <- data.frame(time=test_time, status=test_status, edge=t(test_mats))
        }else if(dim(x)[2] == 1){
          test_x <- x[leftout]
          test_df <- data.frame(time=test_time, status=test_status, edge=t(test_mats), x=test_x)
        }else if(dim(x)[2] > 1){
          test_x <- x[leftout,]
          test_df <- data.frame(time=test_time, status=test_status, edge=t(test_mats), x=test_x)
        }else{
          stop("x should be a DataFrame")
        }
        
        fit <- coxph(Surv(time,status) ~ ., data=train_df)
        
        # predict TEST sub survival
        
        lp_pred[leftout] <- predict(fit, newdata=test_df)
        status_actual[leftout] <- test_status
        time_actual[leftout] <- test_time
        
      }
      
      
    }
    
    if(edge=="separate"){
      
      return(list(positive_edges=pos_edge_index, 
                  negative_edges=neg_edge_index,
                  positive_predicted_linear_predictor=lp_pred_pos, 
                  negative_predicted_linear_predictor=lp_pred_neg,
                  actual_status=status_actual, actual_time=time_actual))
      
    }else if(edge=="combined"){
      
      return(list(positive_edges=pos_edge_index, 
                  negative_edges=neg_edge_index,
                  predicted_linear_predictor=lp_pred, 
                  actual_status=status_actual, actual_time=time_actual))
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
      
      edge_index_pos <- vector(mode="list", length=k)
      edge_index_neg <- vector(mode="list", length=k)
      edge_index <- vector(mode="list", length=k)
        
      lp_pred_pos <- matrix(0,k,samplesize)
      lp_pred_neg <- matrix(0,k,samplesize)
      lp_pred <- matrix(0,k,samplesize)
      status_actual <- matrix(0,k,samplesize)
      time_actual <- matrix(0,k,samplesize)

      for (fold in 1:k){
          
        test_index <- final_index[fold,]
        train_index <- as.vector(final_index[-fold,])
          
        # divide into TRAIN set and TEST set
          
        train_mats <- all_edges[train_index,]
        train_time <- time[train_index]
        train_status <- status[train_index]
        
        test_mats <- all_edges[test_index,]
        test_time <- time[test_index]
        test_status <- status[test_index]
        
        for(j in 1:M){
          ms_edge <- data.frame(time=train_time, status=train_status, edge=train_mats[,j])
          fit_cox <- coxph(Surv(time,status) ~ ., data=ms_edge)
          coef[j] <- summary(fit_cox)$coefficients[1]
          pval[j] <- summary(fit_cox)$coefficients[5]
        }
        
        edge_index_pos[[fold]] <- which(coef > 0 & pval<=thresh)
        edge_index_neg[[fold]] <- which(coef < 0 & pval<=thresh)

        if(edge=="separate"){
          
          train_mats_pos <- rowSums(train_mats[,which(coef > 0 & pval<=thresh)])
          train_mats_neg <- rowSums(train_mats[,which(coef < 0 & pval<=thresh)])
          
          test_mats_pos <- rowSums(test_mats[,which(coef > 0 & pval<=thresh)])
          test_mats_neg <- rowSums(test_mats[,which(coef < 0 & pval<=thresh)])
          
          if(missing(x)){
            train_x <- NULL
            train_df_pos <- data.frame(time=train_time, status=train_status, edge=train_mats_pos)
            train_df_neg <- data.frame(time=train_time, status=train_status, edge=train_mats_neg)
            
          }else if(dim(x)[2] == 1){
            train_x <- x[train_index]
            train_df_pos <- data.frame(time=train_time, status=train_status, edge=train_mats_pos, x=train_x)
            train_df_neg <- data.frame(time=train_time, status=train_status, edge=train_mats_neg, x=train_x)
            
          }else if(dim(x)[2] > 1){
            train_x <- x[train_index,]
            train_df_pos <- data.frame(time=train_time, status=train_status, edge=train_mats_pos, x=train_x)
            train_df_neg <- data.frame(time=train_time, status=train_status, edge=train_mats_neg, x=train_x)
            
          }else{
            stop("x should be a DataFrame")
          }
          
          if(missing(x)){
            test_x <- NULL
            test_df_pos <- data.frame(time=test_time, status=test_status, edge=test_mats_pos)
            test_df_neg <- data.frame(time=test_time, status=test_status, edge=test_mats_neg)
            
          }else if(dim(x)[2] == 1){
            test_x <- x[test_index]
            test_df_pos <- data.frame(time=test_time, status=test_status, edge=test_mats_pos, x=test_x)
            test_df_neg <- data.frame(time=test_time, status=test_status, edge=test_mats_neg, x=test_x)
            
          }else if(dim(x)[2] > 1){
            test_x <- x[test_index,]
            test_df_pos <- data.frame(time=test_time, status=test_status, edge=test_mats_pos, x=test_x)
            test_df_neg <- data.frame(time=test_time, status=test_status, edge=test_mats_neg, x=test_x)
            
          }else{
            stop("x should be a DataFrame")
          }
          
          fit_pos <- coxph(Surv(time,status) ~ .,data=train_df_pos)
          fit_neg <- coxph(Surv(time,status) ~ .,data=train_df_neg)
          
          # predict TEST sub survival
          
          lp_pred_pos[fold,] <- predict(fit_pos, newdata=test_df_pos)
          lp_pred_neg[fold,] <- predict(fit_neg, newdata=test_df_neg)
          status_actual[fold,] <- test_status
          time_actual[fold,] <- test_time
          
        }else if(edge=="combined"){
          
          train_mats <- rowSums(train_mats[,which(pval<=thresh)])
          
          test_mats <- rowSums(test_mats[,which(pval<=thresh)])
          
          if(missing(x)){
            train_x <- NULL
            train_df <- data.frame(time=train_time, status=train_status, edge=train_mats)
          }else if(dim(x)[2] == 1){
            train_x <- x[train_index]
            train_df <- data.frame(time=train_time, status=train_status, edge=train_mats, x=train_x)
          }else if(dim(x)[2] > 1){
            train_x <- x[train_index,]
            train_df <- data.frame(time=train_time, status=train_status, edge=train_mats, x=train_x)
          }else{
            stop("x should be a DataFrame")
          }
          
          if(missing(x)){
            test_x <- NULL
            test_df <- data.frame(time=test_time, status=test_status, edge=test_mats)
          }else if(dim(x)[2] == 1){
            test_x <- x[test_index]
            test_df <- data.frame(time=test_time, status=test_status, edge=test_mats, x=test_x)
          }else if(dim(x)[2] > 1){
            test_x <- x[test_index,]
            test_df <- data.frame(time=test_time, status=test_status, edge=test_mats, x=test_x)
          }else{
            stop("x should be a DataFrame")
          }
          
          fit <- coxph(Surv(time,status) ~ ., data=train_df)
          
          # predict TEST sub survival
          
          lp_pred[fold,] <- predict(fit, newdata=test_df)
          status_actual[fold,] <- test_status
          time_actual[fold,] <- test_time
          
        }
      
          
      }
      
      if(edge=="separate"){
        
        return(list(positive_edges=edge_index_pos, 
                    negative_edges=edge_index_neg,
                    positive_predicted_linear_predictor=lp_pred_pos,
                    negative_predicted_linear_predictor=lp_pred_neg,
                    actual_status=status_actual, 
                    actual_time=time_actual))
        
      }else if(edge=="combined"){
        
        return(list(positive_edges=edge_index_pos, 
                    negative_edges=edge_index_neg,
                    predicted_linear_predictor=lp_pred,
                    actual_status=status_actual, 
                    actual_time=time_actual))
      }
        
    }
      
  
}

}

