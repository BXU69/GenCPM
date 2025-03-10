#' Fit a Linear Model for Continuous Response using Connectome-based Predictive Modeling
#'
#' @import psych
#' @param connectome an array indicating the connectivity between M edges and over N subjects. The dimension should be `M*M*N`.
#' @param behavior a vector containing the behavior measure for all subjects.
#' @param x a data frame containing the non-image variables in the model.
#' @param cv a character indicating the method of cross-validation. The default method is "leave-one-out" cross-validation.
#' @param k a parameter used to set the number of folds for k-fold cross-validation.
#' @param thresh the value of the threshold for selecting significantly related edges. The default value is .01.
#' @param edge a character indicating the model is fitted with either positive and negative edges respectively or combined edges together. The default is "separate".
#' @param seed the value used to set seed for random sampling in the process of cross-validation. The default value is 1220.
#' @return A list contains positive edges, negative edges, Pearson correlation coefficient, p-value, predicted behavior, and actual behavior.
#' @export

linear.GenCPM <- function(connectome, behavior, x=NULL,
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

      # leave out subjects from matrices and behavior

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

      cpm <- train.GenCPM(train_array,train_behav,train_x,fit="linear",thresh=thresh,edge=edge)
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
        stop("edges must be separate or combined")
      }

      # run model on TEST subs

      test_array <- connectome[,,leftout]
      test_behav <- behavior[leftout]

      test_sumpos <- sum(test_array*pos_edges, na.rm = TRUE)/2
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
          behav_pred_pos[leftout] <- predict(fit_pos, newdata=test_pos)
        }else{
          behav_pred_pos[leftout] <- NA
        }

        if((NA %in% fit_neg) == FALSE){
        behav_pred_neg[leftout] <- predict(fit_neg, newdata=test_neg)
        }else{
        behav_pred_neg[leftout] <- NA
        }

      }else if(edge == "combined"){

        if((NA %in% fit_combined) == FALSE){
          behav_pred[leftout] <- predict(fit_combined, newdata=test_combined)
        }
        else{
          behav_pred[leftout] <- NA
        }
        }

        behav_actual[leftout] <- test_behav

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

      rand_index <- seq(1,N,by=1)
      randinds <- sample(rand_index,replace = FALSE) # shuffle random index

      samplesize <- floor(N/k)

      behav_pred_pos <- matrix(0,k,samplesize)
      behav_pred_neg <- matrix(0,k,samplesize)
      behav_pred <- matrix(0,k,samplesize)

      behav_actual <- matrix(0,k,samplesize)

      r_mat <- vector(mode='list', length=k)
      p_mat <- vector(mode='list', length=k)
      positive_edges <- vector(mode='list', length=k)
      negative_edges <- vector(mode='list', length=k)


      for (fold in 1:k){

        s <- 1+(fold-1)*samplesize
        f <- fold*samplesize

        # select one fold as TEST set

        test_index <- randinds[s:f]
        train_index <- randinds[! randinds %in% test_index]

        # divide TRAIN and TEST set

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

        test_array <- connectome[,,test_index] # "samplesize"'s matrices
        test_behav <- behavior[test_index]

        # generate a matrix to save actual behaviors
        behav_actual[fold,] <- test_behav

        # train GenCPM model

        cpm <- train.GenCPM(train_array,train_behav,train_x,fit="linear",thresh=thresh,edge=edge)
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
          stop("edge must be separate or combined")
        }

        # calculate the sum of TEST subs

        test_sumpos <- rep(0,samplesize)
        test_sumneg <- rep(0,samplesize)

        for (j in 1:samplesize){
          test_sumpos[j] <- sum(test_array[,,j]*pos_edges, na.rm = TRUE)/2
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
            behav_pred_pos[fold,] <- predict(fit_pos, newdata=test_pos)  # test_pos%*%as.vector(coef(fit_pos))
          }
          else{
            behav_pred_pos[fold,] <- NA
          }

          if((NA %in% fit_neg) == FALSE){
            behav_pred_neg[fold,] <- predict(fit_neg, newdata=test_neg)
          }
          else{
            behav_pred_neg[fold,] <- NA
          }

        }else if(edge == "combined"){

          if((NA %in% fit_combined) == FALSE){
            behav_pred[fold,] <- predict(fit_combined, newdata=test_combined)
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
