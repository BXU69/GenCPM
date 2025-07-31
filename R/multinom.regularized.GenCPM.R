#' Penalized Version of Multinomial Logistic Model Prediction
#'
#' @import psych
#' @import glmnet
#' @param connectome an array indicating the connectivity between M edges and over N subjects for model fitting. The dimension should be `M*M*N`.
#' @param behavior a vector containing the behavior measure for all subjects for model fitting.
#' @param x a data frame containing the non-image variables for model fitting.
#' @param external.connectome an external array indicating the connectivity for prediction.
#' @param external.x an external data frame containing the non-image variables for prediction.
#' @param cv a character indicating the method of cross-validation. The default method is "leave-one-out" cross-validation.
#' @param k a parameter used to set the number of folds for k-fold cross-validation.
#' @param correlation the method for finding the correlation between edge and behavior. The default is "pearson". Alternative approaches are "spearman" and "kendall".
#' @param thresh the value of the threshold for selecting significantly related edges. The default value is .01.
#' @param edge a character indicating the model is fitted with either positive and negative edges respectively or combined edges together. The default is "separate".
#' @param type type of penalty. The default is lasso.
#' @param lambda the value of penalty of LASSO regression.
#' @param alpha the alpha for elastic net penalty.
#' @param seed the value used to set seed for random sampling in the process of cross-validation. The default value is 1220.
#' @return A list contains positive edges, negative edges, predicted behavior, actual behavior, and the optimal lambda.
#' @export

multinom.regularized.GenCPM <- function(connectome, behavior, x,
                                        external.connectome = NULL, external.x = NULL,
                                     cv="leave-one-out", k = dim(connectome)[3],
                                     correlation = "pearson", thresh = .01,
                                     edge = "separate", type="lasso",
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
    stop("type should only be lasso, ridge or elastic net")
  }

  # parse the dimensions of data
  N <- dim(connectome)[3]
  M1 <- dim(connectome)[1]
  M <- M1*(M1-1)/2

  c <- length(unique(behavior))

  # convert connectome to edge matrix

  edges <- matrix(0, nrow=M,ncol=N) # edges is M*N matrix

  for (i in 1:N){
    edges[,i] <- t(connectome[,,i])[lower.tri(connectome[,,i],diag=FALSE)]
  }

  if(missing(external.connectome)){
    stop("Connectome are required for prediction.")
  }else if(missing(external.connectome) & missing(external.x)){
  # leave-one-out

  if((cv == "leave-one-out" && k == dim(connectome)[3]) | (cv == "k-fold" && k == dim(connectome)[3])){

    behav_pred <- rep(0,N)
    behav_pred_pos <- rep(0,N)
    behav_pred_neg <- rep(0,N)

    behav_actual <- rep(0,N)

    lambda_total <- rep(0,N)
    lambda_total_pos <- rep(0,N)
    lambda_total_neg <- rep(0,N)

    selected_edges_pos <- vector(mode='list', length=N)
    selected_edges_neg <- vector(mode='list', length=N)

    for(leftout in 1:N){

      # divide into TRAIN set and TEST set

      train_mats <- edges[,-leftout]
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

      test_mats <- edges[,leftout]
      test_behav <- behavior[leftout]

      if(missing(x)){
        test_x <- NULL
      }else if(dim(x)[2] == 1){
        test_x <- x[leftout]
      }else if(dim(x)[2] > 1){
        test_x <- x[leftout,]
      }else{
        stop("x should be a DataFrame")
      }

      # univariate edge selection
      # correlate all edges with behavior

      corr <- corr.test(x = t(train_mats), y = train_behav, method = correlation,
                        adjust = "none", ci=F)

      edge_r <- corr$r
      edge_p <- corr$p

      # set threshold and define masks

      edge_index <- which(edge_p <= thresh)

      edge_index_pos <- which(edge_r > 0 & edge_p <= thresh)
      edge_index_neg <- which(edge_r < 0 & edge_p <= thresh)

      train_covars <- data.matrix(cbind(edge=t(train_mats[edge_index,]), x=train_x))
      test_covars <- data.matrix(cbind(edge=t(test_mats[edge_index]), x=test_x))

      train_covars_pos <- data.matrix(cbind(edge=t(train_mats[edge_index_pos,]), x=train_x))
      test_covars_pos <- data.matrix(cbind(edge=t(test_mats[edge_index_pos]), x=test_x))

      train_covars_neg <- data.matrix(cbind(edge=t(train_mats[edge_index_neg,]), x=train_x))
      test_covars_neg <- data.matrix(cbind(edge=t(test_mats[edge_index_neg]), x=test_x))

      nedge <- length(edge_index)
      nedge_pos <- length(edge_index_pos)
      nedge_neg <- length(edge_index_neg)

      # lambda selection

      if(edge=="separate"){

        if(missing(lambda) | length(lambda) < 1){

          fit_lambda_pos <- cv.glmnet(train_covars_pos, train_behav,
                                      family="multinomial", alpha=type_alpha,
                                      type.measure="class",
                                      penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
          opt_lambda_pos <- fit_lambda_pos$lambda.min

          fit_lambda_neg <- cv.glmnet(train_covars_neg, train_behav,
                                      family="multinomial", alpha=type_alpha,
                                      type.measure="class",
                                      penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
          opt_lambda_neg <- fit_lambda_neg$lambda.min

        }else if(length(lambda) > 1){

          fit_lambda_pos <- cv.glmnet(train_covars_pos, train_behav,
                                      family="multinomial", lambda = lambda,
                                      alpha=type_alpha, type.measure="class",
                                      penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
          opt_lambda_pos <- fit_lambda_pos$lambda.min

          fit_lambda_neg <- cv.glmnet(train_covars_neg, train_behav,
                                      family="multinomial", lambda = lambda,
                                      alpha=type_alpha, type.measure="class",
                                      penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
          opt_lambda_neg <- fit_lambda_neg$lambda.min


        }else{

          opt_lambda_pos <- lambda
          opt_lambda_neg <- lambda

        }

        lambda_total_pos[leftout] <- opt_lambda_pos
        lambda_total_neg[leftout] <- opt_lambda_neg

        # train model with the optimal lambda

        fit_pos <- glmnet(train_covars_pos, train_behav, family="multinomial",
                          lambda = opt_lambda_pos, alpha=type_alpha,
                          penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
        fit_neg <- glmnet(train_covars_neg, train_behav, family="multinomial",
                          lambda = opt_lambda_neg, alpha=type_alpha,
                          penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))

        # edge selection by glmnet

        coef_pos <- matrix(NA,nrow=nedge_pos,ncol=c)
        coef_neg <- matrix(NA,nrow=nedge_neg,ncol=c)

        index0_pos <- vector(mode="list",length=c)
        index0_neg <- vector(mode="list",length=c)

        for (p in 1:c){
          coef_pos[,p] <- data.matrix(coef(fit_pos)[[p]])[c(2:(nedge_pos+1))]
          index0_pos[[p]] <- which(coef_pos[,p]!=0)
          coef_neg[,p] <- data.matrix(coef(fit_neg)[[p]])[c(2:(nedge_neg+1))]
          index0_neg[[p]] <- which(coef_neg[,p]!=0)
        }

        index1_pos <- sort(unique(unlist(index0_pos)))
        index1_neg <- sort(unique(unlist(index0_neg)))

        selected_edges_pos[[leftout]] <- edge_index_pos[index1_pos]
        selected_edges_neg[[leftout]] <- edge_index_neg[index1_neg]

        # predict TEST subs with the best lambda param

        behav_pred_pos[leftout] <- predict(fit_pos, newx=test_covars_pos, s = opt_lambda_pos, type="class")
        behav_pred_neg[leftout] <- predict(fit_neg, newx=test_covars_neg, s = opt_lambda_neg, type="class")

        behav_actual[leftout] <- test_behav


      }else if(edge=="combined"){

        if(missing(lambda) | length(lambda) < 1){

          fit_lambda <- cv.glmnet(train_covars, train_behav,
                                  family="multinomial", alpha=type_alpha,
                                  type.measure="class",
                                  penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
          opt_lambda <- fit_lambda$lambda.min

        }else if(length(lambda) > 1){

          fit_lambda <- cv.glmnet(train_covars, train_behav,
                                  family="multinomial", lambda = lambda,
                                  alpha=type_alpha, type.measure="class",
                                  penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
          opt_lambda <- fit_lambda$lambda.min

        }else{

          opt_lambda <- lambda

        }

        lambda_total[leftout] <- opt_lambda

        # train model with the optimal lambda

        fit <- glmnet(train_covars, train_behav, family="multinomial",
                      lambda = opt_lambda, alpha=type_alpha,
                      penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))

        # edge selection by glmnet

        coef <- matrix(NA,nrow=nedge,ncol=c)
        index0_pos <- vector(mode="list",length=c)
        index0_neg <- vector(mode="list",length=c)

        for (p in 1:c){
          coef[,p] <- data.matrix(coef(fit)[[p]])[c(2:(nedge+1))]
          index0_pos[[p]] <- which(coef[,p]>0)
          index0_neg[[p]] <- which(coef[,p]<0)
        }

        index1_pos <- sort(unique(unlist(index0_pos)))
        index1_neg <- sort(unique(unlist(index0_neg)))

        selected_edges_pos[[leftout]] <- edge_index[index1_pos]
        selected_edges_neg[[leftout]] <- edge_index[index1_neg]

        # predict TEST sub with the best lambda param

        behav_pred[leftout] <- predict(fit, newx=test_covars, s = opt_lambda, type="class")
        behav_actual[leftout] <- test_behav

      }

    }

    if(edge=="separate"){

      return(list(positive_edges=selected_edges_pos,
                  negative_edges=selected_edges_neg,
                  positive_predicted_behavior=behav_pred_pos,
                  negative_predicted_behavior=behav_pred_neg,
                  actual_behavior=behav_actual,
                  positive_lambda_total=lambda_total_pos,
                  negative_lambda_total=lambda_total_neg))

    }
    else if(edge=="combined"){

      return(list(positive_edges=selected_edges_pos,
                  negative_edges=selected_edges_neg,
                  predicted_behavior=behav_pred,
                  actual_behavior=behav_actual,
                  lambda_total=lambda_total))

    }

  }

  if(cv == "k-fold" | ((cv == "leave-one-out") && (k != dim(connectome)[3]))){

    if(k-floor(k) != 0){
      stop("k should be an integer.")
    }
    else{

      set.seed(seed)

      labels <- unique(behavior)
      ratio <- rep(NA,length(labels))

      for (i in 1:length(labels)){
        ratio[i] <- length(which(behavior==labels[i]))/N
      }

      n_label <- floor(N/k *ratio)
      samplesize <- sum(n_label)

      index_list <- vector(mode="list",length=length(labels))
      randinds <- vector(mode="list",length=length(labels))
      final_index <- matrix(NA,nrow=k,ncol=samplesize)

      for (i in 1:length(labels)){
        index_list[[i]] <- which(behavior==labels[[i]])
      }

      for (fold in 1:k){

        for(i in 1:length(labels)){
          randinds[[i]] <- sample(index_list[[i]], size=n_label[i], replace=FALSE)
        }

        final_index[fold,] <- sample(unlist(randinds), size=samplesize, replace=FALSE)

        for (i in 1:length(labels)){
          index_list[[i]] <- setdiff(index_list[[i]],randinds[[i]])
        }

      }

      behav_pred <- matrix(0,k,samplesize)
      behav_pred_pos <- matrix(0,k,samplesize)
      behav_pred_neg <- matrix(0,k,samplesize)

      behav_actual <- matrix(0,k,samplesize)

      lambda_total <- rep(0,k)
      lambda_total_pos <- rep(0,k)
      lambda_total_neg <- rep(0,k)

      selected_edges_pos <- vector(mode='list', length=k)
      selected_edges_neg <- vector(mode='list', length=k)

      for (fold in 1:k){

        test_index <- final_index[fold,]
        train_index <- as.vector(final_index[-fold,])

        # divide into TRAIN set and TEST set

        train_mats <- edges[,train_index]
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

        test_mats <- edges[,test_index]
        test_behav <- behavior[test_index]

        if(missing(x)){
          test_x <- NULL
        }else if(dim(x)[2] == 1){
          test_x <- x[test_index]
        }else if(dim(x)[2] > 1){
          test_x <- x[test_index,]
        }else{
          stop("x should be a DataFrame")
        }

        # univariate edge selection
        # correlate all edges with behavior

        corr <- corr.test(x = t(train_mats), y = train_behav, method = correlation,
                          adjust = "none", ci=F)

        edge_r <- corr$r
        edge_p <- corr$p

        # set threshold and define masks

        edge_index <- which(edge_p <= thresh)

        edge_index_pos <- which(edge_r > 0 & edge_p <= thresh)
        edge_index_neg <- which(edge_r < 0 & edge_p <= thresh)

        train_covars <- data.matrix(cbind(edge=t(train_mats[edge_index,]), x=train_x))
        test_covars <- data.matrix(cbind(edge=t(test_mats[edge_index]), x=test_x))

        train_covars_pos <- data.matrix(cbind(edge=t(train_mats[edge_index_pos,]), x=train_x))
        test_covars_pos <- data.matrix(cbind(edge=t(test_mats[edge_index_pos]), x=test_x))

        train_covars_neg <- data.matrix(cbind(edge=t(train_mats[edge_index_neg,]), x=train_x))
        test_covars_neg <- data.matrix(cbind(edge=t(test_mats[edge_index_neg]), x=test_x))

        nedge <- length(edge_index)
        nedge_pos <- length(edge_index_pos)
        nedge_neg <- length(edge_index_neg)

        if(edge=="separate"){

          if(missing(lambda) | length(lambda) < 1){

            fit_lambda_pos <- cv.glmnet(train_covars_pos, train_behav,
                                        family="multinomial", alpha=type_alpha,
                                        type.measure="class",
                                        penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
            opt_lambda_pos <- fit_lambda_pos$lambda.1se

            fit_lambda_neg <- cv.glmnet(train_covars_neg, train_behav,
                                        family="multinomial", alpha=type_alpha,
                                        type.measure="class",
                                        penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
            opt_lambda_neg <- fit_lambda_neg$lambda.1se

          }else if(length(lambda) > 1){

            fit_lambda_pos <- cv.glmnet(train_covars_pos, train_behav,
                                        family="multinomial", lambda = lambda,
                                        alpha=type_alpha,
                                        type.measure="class",
                                        penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
            opt_lambda_pos <- fit_lambda_pos$lambda.1se

            fit_lambda_neg <- cv.glmnet(train_covars_neg, train_behav,
                                        family="multinomial", lambda = lambda,
                                        alpha=type_alpha,
                                        type.measure="class",
                                        penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
            opt_lambda_neg <- fit_lambda_neg$lambda.1se

          }else{

            opt_lambda_pos <- lambda
            opt_lambda_neg <- lambda

          }

          lambda_total_pos[fold] <- opt_lambda_pos
          lambda_total_neg[fold] <- opt_lambda_neg

          # train model with the optimal lambda

          fit_pos <- glmnet(train_covars_pos, train_behav, family="multinomial",
                            lambda = opt_lambda_pos, alpha=type_alpha,
                            penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
          fit_neg <- glmnet(train_covars_neg, train_behav, family="multinomial",
                            lambda = opt_lambda_neg, alpha=type_alpha,
                            penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))

          # edge selection by glmnet

          coef_pos <- matrix(NA,nrow=nedge_pos,ncol=c)
          coef_neg <- matrix(NA,nrow=nedge_neg,ncol=c)

          index0_pos <- vector(mode="list",length=c)
          index0_neg <- vector(mode="list",length=c)

          for (p in 1:c){
            coef_pos[,p] <- data.matrix(coef(fit_pos)[[p]])[c(2:(nedge_pos+1))]
            index0_pos[[p]] <- which(coef_pos[,p]!=0)
            coef_neg[,p] <- data.matrix(coef(fit_neg)[[p]])[c(2:(nedge_neg+1))]
            index0_neg[[p]] <- which(coef_neg[,p]!=0)
          }

          index1_pos <- sort(unique(unlist(index0_pos)))
          index1_neg <- sort(unique(unlist(index0_neg)))

          selected_edges_pos[[fold]] <- edge_index_pos[index1_pos]
          selected_edges_neg[[fold]] <- edge_index_neg[index1_neg]

          # predict TEST subs with the best lambda param

          behav_pred_pos[fold,] <- predict(fit_pos, newx=test_covars_pos, s = opt_lambda_pos, type="class")
          behav_pred_neg[fold,] <- predict(fit_neg, newx=test_covars_neg, s = opt_lambda_neg, type="class")

          behav_actual[fold,] <- test_behav


        }else if(edge=="combined"){

          if(missing(lambda) | length(lambda) < 1){

            fit_lambda <- cv.glmnet(train_covars, train_behav,
                                    family="multinomial", alpha=type_alpha,
                                    type.measure="class",
                                    penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
            opt_lambda <- fit_lambda$lambda.1se

          }else if(length(lambda) > 1){

            fit_lambda <- cv.glmnet(train_covars, train_behav,
                                    family="multinomial", lambda = lambda,
                                    alpha=type_alpha,
                                    type.measure="class",
                                    penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
            opt_lambda <- fit_lambda$lambda.1se

          }else{

            opt_lambda <- lambda

          }

          lambda_total[fold] <- opt_lambda

          # train model with the optimal lambda

          fit <- glmnet(train_covars, train_behav, family="multinomial",
                        lambda = opt_lambda, alpha=type_alpha,
                        penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))

          # edge selection by glmnet

          coef <- matrix(NA,nrow=nedge,ncol=c)
          index0_pos <- vector(mode="list",length=c)
          index0_neg <- vector(mode="list",length=c)

          for (p in 1:c){
            coef[,p] <- data.matrix(coef(fit)[[p]])[c(2:(nedge+1))]
            index0_pos[[p]] <- which(coef[,p]>0)
            index0_neg[[p]] <- which(coef[,p]<0)
          }

          index1_pos <- sort(unique(unlist(index0_pos)))
          index1_neg <- sort(unique(unlist(index0_neg)))

          selected_edges_pos[[fold]] <- edge_index[index1_pos]
          selected_edges_neg[[fold]] <- edge_index[index1_neg]

          # predict TEST subs with the best lambda param

          behav_pred[fold,] <- predict(fit, newx=test_covars, s = opt_lambda, type="class")
          behav_actual[fold,] <- test_behav

        }

      }

      if(edge=="separate"){

        return(list(positive_edges=selected_edges_pos,
                    negative_edges=selected_edges_neg,
                    positive_predicted_behavior=behav_pred_pos,
                    negative_predicted_behavior=behav_pred_neg,
                    actual_behavior=behav_actual,
                    positive_lambda_total=lambda_total_pos,
                    negative_lambda_total=lambda_total_neg))

      }else if(edge=="combined"){

        return(list(positive_edges=selected_edges_pos,
                    negative_edges=selected_edges_neg,
                    predicted_behavior=behav_pred,
                    actual_behavior=behav_actual,
                    lambda_total=lambda_total))

      }

    }

  }


  }else{
    
    # divide into TRAIN set and TEST set
    
    train_mats <- edges
    train_behav <- behavior
    
    if(missing(x)){
      train_x <- NULL
    }else {
      train_x <- x
    }
    
    t <- dim(external.connectome)[3]
    
    external_edges <- matrix(0, nrow=M,ncol=t)
    
    for (i in 1:t){
      external_edges[,i] <- t(connectome[,,i])[lower.tri(connectome[,,i],diag=FALSE)]
    }
    
    test_mats <- external_edges
    
    if(missing(x)){
      test_x <- NULL
    }else {
      test_x <- external.x
    }
    
    # univariate edge selection
    
    corr <- corr.test(x = t(train_mats), y = train_behav, method = correlation,
                      adjust = "none", ci=F)
    
    edge_r <- as.vector(corr$r)
    edge_p <- as.vector(corr$p)
    
    # set threshold and define masks
    
    edge_index <- which(edge_p <= thresh)
    
    edge_index_pos <- which(edge_r > 0 & edge_p <= thresh)
    edge_index_neg <- which(edge_r < 0 & edge_p <= thresh)
    
    train_covars <- data.matrix(cbind(edge=t(train_mats[edge_index,]), x=train_x))
    test_covars <- data.matrix(cbind(edge=t(test_mats[edge_index]), x=test_x))
    
    train_covars_pos <- data.matrix(cbind(edge=t(train_mats[edge_index_pos,]), x=train_x))
    test_covars_pos <- data.matrix(cbind(edge=t(test_mats[edge_index_pos]), x=test_x))
    
    train_covars_neg <- data.matrix(cbind(edge=t(train_mats[edge_index_neg,]), x=train_x))
    test_covars_neg <- data.matrix(cbind(edge=t(test_mats[edge_index_neg]), x=test_x))
    
    nedge <- length(edge_index)
    nedge_pos <- length(edge_index_pos)
    nedge_neg <- length(edge_index_neg)
    
    if(edge=="separate"){
      
      if(missing(lambda) | length(lambda) < 1){
        
        fit_lambda_pos <- cv.glmnet(train_covars_pos, train_behav,
                                    family = "multinomial",
                                    alpha=type_alpha, type.measure="class",
                                    penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
        opt_lambda_pos <- fit_lambda_pos$lambda.1se
        
        fit_lambda_neg <- cv.glmnet(train_covars_neg, train_behav,
                                    family = "multinomial",
                                    alpha=type_alpha, type.measure="class",
                                    penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
        opt_lambda_neg <- fit_lambda_neg$lambda.1se
        
      }else if(length(lambda) > 1){
        
        fit_lambda_pos <- cv.glmnet(train_covars_pos, train_behav,
                                    family = "multinomial",
                                    lambda = lambda, alpha=type_alpha,
                                    type.measure="class",
                                    penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
        opt_lambda_pos <- fit_lambda_pos$lambda.1se
        
        fit_lambda_neg <- cv.glmnet(train_covars_neg, train_behav,
                                    family = "multinomial",
                                    lambda = lambda, alpha=type_alpha,
                                    type.measure="class",
                                    penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
        opt_lambda_neg <- fit_lambda_neg$lambda.1se
        
      }else{
        
        opt_lambda_pos <- lambda
        opt_lambda_neg <- lambda
        
      }
      
      lambda_total_pos <- opt_lambda_pos
      lambda_total_neg <- opt_lambda_neg
      
      # train model with the optimal lambda
      
      fit_pos <- glmnet(train_covars_pos, train_behav,
                        family = "multinomial",
                        lambda = opt_lambda_pos, alpha=type_alpha,
                        penalty.factor=c(rep(1, nedge_pos),rep(0, ncol(train_covars_pos)-nedge_pos)))
      fit_neg <- glmnet(train_covars_neg, train_behav,
                        family = "multinomial",
                        lambda = opt_lambda_neg, alpha=type_alpha,
                        penalty.factor=c(rep(1, nedge_neg),rep(0, ncol(train_covars_neg)-nedge_neg)))
      
      # edge selection by glmnet
      
      coef_pos <- matrix(NA,nrow=nedge_pos,ncol=c)
      coef_neg <- matrix(NA,nrow=nedge_neg,ncol=c)
      
      index0_pos <- vector(mode="list",length=c)
      index0_neg <- vector(mode="list",length=c)
      
      for (p in 1:c){
        coef_pos[,p] <- data.matrix(coef(fit_pos)[[p]])[c(2:(nedge_pos+1))]
        index0_pos[[p]] <- which(coef_pos[,p]!=0)
        coef_neg[,p] <- data.matrix(coef(fit_neg)[[p]])[c(2:(nedge_neg+1))]
        index0_neg[[p]] <- which(coef_neg[,p]!=0)
      }
      
      index1_pos <- sort(unique(unlist(index0_pos)))
      index1_neg <- sort(unique(unlist(index0_neg)))
      
      selected_edges_pos <- edge_index_pos[index1_pos]
      selected_edges_neg <- edge_index_neg[index1_neg]
      
      # predict TEST subs with the best lambda param
      
      behav_pred_pos <- predict(fit_pos, newx=test_covars_pos, type = "class")
      behav_pred_neg <- predict(fit_neg, newx=test_covars_neg, type = "class")
      
      
    }else if(edge=="combined"){
      
      if(missing(lambda) | length(lambda) < 1){
        
        fit_lambda <- cv.glmnet(train_covars, train_behav,
                                family = "multinomial",
                                alpha=type_alpha, type.measure="class",
                                penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
        opt_lambda <- fit_lambda$lambda.1se
        
      }else if(length(lambda) > 1){
        
        fit_lambda <- cv.glmnet(train_covars, train_behav,
                                family = "multinomial",
                                lambda = lambda, alpha=type_alpha,
                                type.measure="class",
                                penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
        opt_lambda <- fit_lambda$lambda.1se
        
      }else{
        
        opt_lambda <- lambda
        
      }
      
      lambda_total <- opt_lambda
      
      # train model with the optimal lambda
      
      fit <- glmnet(train_covars, train_behav, 
                    family = "multinomial", lambda = opt_lambda,
                    alpha=type_alpha,
                    penalty.factor=c(rep(1, nedge),rep(0, ncol(train_covars)-nedge)))
      
      # edge selection by glmnet
      
      coef <- matrix(NA,nrow=nedge,ncol=c)
      index0_pos <- vector(mode="list",length=c)
      index0_neg <- vector(mode="list",length=c)
      
      for (p in 1:c){
        coef[,p] <- data.matrix(coef(fit)[[p]])[c(2:(nedge+1))]
        index0_pos[[p]] <- which(coef[,p]>0)
        index0_neg[[p]] <- which(coef[,p]<0)
      }
      
      index1_pos <- sort(unique(unlist(index0_pos)))
      index1_neg <- sort(unique(unlist(index0_neg)))
      
      selected_edges_pos <- edge_index[index1_pos]
      selected_edges_neg <- edge_index[index1_neg]
      
      
      # predict TEST subs with the best lambda param
      
      behav_pred <- predict(fit, newx=test_covars, type = "class")
      
    }
    
    
    if(edge=="separate"){
      
      return(list(positive_edges=selected_edges_pos,
                  negative_edges=selected_edges_neg,
                  positive_model=fit_pos,
                  negative_model=fit_neg,
                  positive_predicted_behavior=behav_pred_pos,
                  negative_predicted_behavior=behav_pred_neg,
                  positive_lambda=lambda_total_pos,
                  negative_lambda=lambda_total_neg))
      
    }else if(edge=="combined"){
      
      return(list(positive_edges=selected_edges_pos,
                  negative_edges=selected_edges_neg,
                  combined_model=fit,
                  predicted_behavior=behav_pred,
                  lambda=lambda_total))
      
    }
    
  }
  
  
}








