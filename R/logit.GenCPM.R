#' Fit a Logistic Model for Binary Response using Connectome-based Predictive Modeling
#' 
#' @import psych
#' @param connectome a array indicating the connectivity between M edges and over N subjects. The dimension should be `M*M*N`.
#' @param behavior a vector containing the behavior measure for all subjects.
#' @param x a data frame containing the non-image variables in the model.
#' @param cv a character indicating the method of cross-validation. The default method is "leave-one-out" cross-validation.
#' @param k a parameter used to set the number of folds for k-fold cross-validation.
#' @param thresh the value of the threshold for selecting significantly related edges. The default value is .01.
#' @param edge a character indicating the model is fitted with either positive and negative edges respectively or combined edges together. The default is "separate".
#' @param seed the value used to set seed for random sampling in the process of cross-validation. The default value is 1220.
#' @return A list contains positive edges, negative edges, Pearson correlation coefficient, p-value, predicted behavior, and actual behavior.
#' @export

logit.GenCPM <- function(connectome, behavior, x=NULL, 
                       cv="leave-one-out", k = dim(connectome)[3], 
                       thresh = .01, edge = "separate", seed = 1220){
  
  M <- dim(connectome)[1]
  N <- dim(connectome)[3]
  
  if(cv != "leave-one-out" && cv != "k-fold"){
    stop("The input of cv can only be leave-one-out or k-fold.")
  }
  
  if((cv == "leave-one-out" && k == dim(connectome)[3]) | (cv == "k-fold" && k == dim(connectome)[3])){
    
    behav_pred_pos <- rep(0, N)
    behav_pred_neg <- rep(0, N)
    behav_pred <- rep(0,N)
    behav_actual <- rep(0,N)
    r_mat <- vector(mode='list', length=N)
    p_mat <- vector(mode='list', length=N)
    positive_edges <- vector(mode='list', length=N)
    negative_edges <- vector(mode='list', length=N)
    
    
    for (leftout in 1:N){
      
      train_array <- connectome[,,-leftout]
      train_behav <- behavior[-leftout]
      
      if(missing(x)){
        train_x <- NULL
      }else if(dim(x)[2] == 1){
        train_x <- x[-leftout]
      }else if(dim(x)[2] > 1){
        train_x <- x[-leftout,]
      }else{
        stop("x should be a DataFrame")
      }
      
      # train GenCPM model
      
      cpm <- train.GenCPM(train_array,train_behav,train_x,fit="logistic",thresh=thresh,edge=edge)
      r_mat[[leftout]] <- cpm$r_mat
      p_mat[[leftout]] <- cpm$p_mat
      pos_edges <- cpm$positive_edges_matrix
      neg_edges <- cpm$negative_edges_matrix
      positive_edges[[leftout]] <- cpm$positive_edges
      negative_edges[[leftout]] <- cpm$negative_edges
      
      if(edge == "separate"){
        fit_pos <- cpm$positive_model
        fit_neg <- cpm$negative_model
      }else if(edge == "combined"){
        fit_combined <- cpm$combined_model
      }else{
        stop("model can only be fitted using either separate or combined edges")
      }
      
      # run model on TEST subs
      
      test_array <- connectome[,,leftout]
      test_behav <- behavior[leftout]
      test_sumpos <- sum(test_array*pos_edges, na.rm = TRUE)/2 # sum(test_array%*%pos_edges[,], na.rm = TRUE)/2
      test_sumneg <- sum(test_array*neg_edges, na.rm = TRUE)/2
      
      if(missing(x)){
        test_x <- NULL
        test_pos <- data.frame(sumpos=test_sumpos)
        test_neg <- data.frame(sumneg=test_sumneg)
        test_combined <- data.frame(sumpos=test_sumpos, sumneg=test_sumneg)
      }else if(dim(x)[2] == 1){
        test_x <- x[leftout]
        test_pos <- data.frame(sumpos=test_sumpos, x=test_x)
        test_neg <- data.frame(sumneg=test_sumneg, x=test_x)
        test_combined <- data.frame(sumpos=test_sumpos, sumneg=test_sumneg, x=test_x)
      }else if(dim(x)[2] > 1){
        test_x <- x[leftout,]
        test_pos <- data.frame(sumpos=test_sumpos, x=test_x)
        test_neg <- data.frame(sumneg=test_sumneg, x=test_x)
        test_combined <- data.frame(sumpos=test_sumpos, sumneg=test_sumneg, x=test_x)
        
      }
      
      if(edge=="separate"){
        if((NA %in% fit_pos) == FALSE){
          behav_pred_pos[leftout] <- predict(fit_pos, newdata=test_pos, type="response")
        }
        else{
          behav_pred_pos[leftout] <- NA
        }
        
        if((NA %in% fit_neg) == FALSE){
          behav_pred_neg[leftout] <- predict(fit_neg, newdata=test_neg, type="response")
        }
        else{
          behav_pred_neg[leftout] <- NA
        }

          behav_actual[leftout] <- test_behav
        
      }else if(edge == "combined"){
        
        if((NA %in% fit_combined) == FALSE){
          behav_pred[leftout] <- predict(fit_combined, newdata=test_combined, type="response")
        }
        else{
          behav_pred[leftout] <- NA
        }
          
          behav_actual[leftout] <- test_behav
        }
        
    }
    
    if(edge=="separate"){
      return(list(positive_edges=positive_edges, 
                  negative_edges=negative_edges, 
                  r_mat=r_mat, p_mat=p_mat, 
                  positive_predicted_behavior=behav_pred_pos, 
                  negative_predicted_behavior=behav_pred_neg, 
                  actual_behavior=behav_actual))
    }else if(edge=="combined"){
      return(list(positive_edges=positive_edges, 
                  negative_edges=negative_edges, 
                  r_mat=r_mat, p_mat=p_mat, 
                  predicted_behavior=behav_pred, 
                  actual_behavior=behav_actual))
    }
    
    
  }
  
  if(cv == "k-fold" | ((cv == "leave-one-out") && (k != dim(connectome)[3]))){
    
    if(k-floor(k) != 0){
      stop("k should be an integer.")
    }
    else{
      
      set.seed(seed)
      
      r0 <- length(which(behavior==0))/N
      r1 <- length(which(behavior==1))/N
      
      n0 <- floor(N/k * r0)
      n1 <- floor(N/k * r1)
      
      samplesize <- n0+n1
      
      index_list <- vector(mode="list",length=2)
      final_index <- matrix(NA,nrow=k,ncol=samplesize)
      
      index_list[[1]] <- which(behavior==0)
      index_list[[2]] <- which(behavior==1)
      
      for (fold in 1:k){
        
        randinds0 <- sample(index_list[[1]], size=n0, replace=FALSE)
        randinds1 <- sample(index_list[[2]], size=n1, replace=FALSE)
        
        final_index[fold,] <- sample(c(randinds0, randinds1),size=samplesize)
        
        index_list[[1]] <- setdiff(index_list[[1]],randinds0)
        index_list[[2]] <- setdiff(index_list[[2]],randinds1)
        
      }
      
      behav_pred_pos <- matrix(0,k,samplesize)
      behav_pred_neg <- matrix(0,k,samplesize)
      behav_pred <- matrix(0,k,samplesize)
      
      behav_actual <- matrix(0,k,samplesize)
      
      r_mat <- vector(mode='list', length=k)
      p_mat <- vector(mode='list', length=k)
      positive_edges <- vector(mode='list', length=k)
      negative_edges <- vector(mode='list', length=k)
      
      for (fold in 1:k){
        
        test_index <- final_index[fold,]
        train_index <- as.vector(final_index[-fold,])
        
        # divide TRAIN set and TEST set
        
        train_array <- connectome[,,train_index]
        train_behav <- behavior[train_index]
        
        if(missing(x)){
          train_x <- NULL
        }else if(dim(x)[2] == 1){
          train_x <- x[train_index]
        }else if(dim(x)[2] > 1){
          train_x <- x[train_index,]
        }else{
          stop("x should be a DataFrame")
        }
        
        test_array <- connectome[,,test_index]
        test_behav <- behavior[test_index]

        behav_actual[fold,] <- test_behav 
        
        # train GenCPM model
        
        cpm <- train.GenCPM(train_array,train_behav,train_x,fit="logistic",thresh=thresh,edge=edge)
        r_mat[[fold]] <- cpm$r_mat
        p_mat[[fold]] <- cpm$p_mat
        pos_edges <- cpm$positive_edges_matrix
        neg_edges <- cpm$negative_edges_matrix
        positive_edges[[fold]] <- cpm$positive_edges
        negative_edges[[fold]] <- cpm$negative_edges
        
        if(edge=="separate"){
          
          fit_pos <- cpm$positive_model
          fit_neg <- cpm$negative_model
          
        }else if(edge == "combined"){
          
          fit_combined <- cpm$combined_model
          
        }else{
          stop("model can only be fitted using either separate or combined edges")
        }
        
        # calculate the sum of TEST subs
        
        test_sumpos <- rep(0,samplesize)
        test_sumneg <- rep(0,samplesize)
        
        for (j in 1:samplesize){
          test_sumpos[j] <- sum(test_array[,,j]*pos_edges, na.rm = TRUE)/2 # sum(test_array[,,j]%*%pos_edges[,], na.rm = TRUE)/2
          test_sumneg[j] <- sum(test_array[,,j]*neg_edges, na.rm = TRUE)/2
        }
        
        if(missing(x)){
          test_x <- NULL
          test_pos <- data.frame(sumpos=test_sumpos)
          test_neg <- data.frame(sumneg=test_sumneg)
          test_combined <- data.frame(sumpos=test_sumpos, sumneg=test_sumneg)
        }else if(dim(x)[2] == 1){
          test_x <- x[test_index]
          test_pos <- data.frame(sumpos=test_sumpos, x=test_x)
          test_neg <- data.frame(sumneg=test_sumneg, x=test_x)
          test_combined <- data.frame(sumpos=test_sumpos, sumneg=test_sumneg, x=test_x)
        }else if(dim(x)[2] > 1){
          test_x <- x[test_index,]
          test_pos <- data.frame(sumpos=test_sumpos, x=test_x)
          test_neg <- data.frame(sumneg=test_sumneg, x=test_x)
          test_combined <- data.frame(sumpos=test_sumpos, sumneg=test_sumneg, x=test_x)
        }
        
        # run model on TEST subs
        
        if(edge == "separate"){
          
          if((NA %in% fit_pos) == FALSE){
            behav_pred_pos[fold,] <- predict(fit_pos, newdata=test_pos, type="response")  # test_pos%*%as.vector(coef(fit_pos))
          }
          else{
            behav_pred_pos[fold,] <- NA
          }
          
          if((NA %in% fit_neg) == FALSE){
              behav_pred_neg[fold,] <- predict(fit_neg, newdata=test_neg, type="response")
          }
          else{
            behav_pred_neg[fold,] <- NA
          }
          
        }else if(edge == "combined"){
          
          if((NA %in% fit_combined) == FALSE){
              behav_pred[fold,] <- predict(fit_combined, newdata=test_combined, type="response")  # test_pos%*%as.vector(coef(fit_pos))
          }
          else{
            behav_pred[fold,] <- NA
          }
          
        }
        
        
      }
      
      if(edge=="separate"){
        return(list(positive_edges=positive_edges, 
                    negative_edges=negative_edges, 
                    r_mat=r_mat, p_mat=p_mat, 
                    positive_predicted_behavior=behav_pred_pos, 
                    negative_predicted_behavior=behav_pred_neg, 
                    actual_behavior=behav_actual))
      }else if(edge=="combined"){
        return(list(positive_edges=positive_edges, 
                    negative_edges=negative_edges, 
                    r_mat=r_mat, p_mat=p_mat, 
                    predicted_behavior=behav_pred, 
                    actual_behavior=behav_actual))
      }
      
    }
    
  }
  
}
    
      
      
   
      
      

  
