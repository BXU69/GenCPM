#' Model Training using Connectome-based Predictive Modeling
#'
#' @import psych
#' @import nnet
#' @param train_array an array indicating the connectivity between M edges and over N subjects. The dimension should be `M*M*N`.
#' @param train_behav a vector containing the behavior measure for all subjects.
#' @param train_x a data frame containing the non-image variables in the model.
#' @param fit the method to be used in fitting the model. The default method is "linear".
#' @param thresh the value of the threshold for selecting significantly related edges. The default value is .01.
#' @param edge a character indicating the model is fitted with either positive and negative edges respectively or combined edges together. The default is "separate".
#' @return A list contains positive edges, negative edges, the model based on positive edges and negative edges separately, or the model based on combined edges.
#' @export

train.GenCPM <- function(train_array, train_behav=NULL, train_x = NULL, fit = "linear", thresh=.01, edge = "separate"){

  M <- dim(train_array)[1]
  N <- dim(train_array)[3]

  # reshape array to matrix

  train_mats <- matrix(train_array, M*M, N)

  # correlate all edges with behavior

  r_mat <- rep(0, M*M)
  p_mat <- rep(0, M*M)

  corr <- corr.test(x = t(train_mats), y = train_behav, adjust = "none", ci=F)

  r_mat <- as.vector(corr$r)
  p_mat <- as.vector(corr$p)

  r_mat <- matrix(r_mat, M, M)
  p_mat <- matrix(p_mat, M, M)

  # set threshold and define masks

  pos_edges <- matrix(rep(0, M*M), M, M)
  neg_edges <- matrix(rep(0, M*M), M, M)

  pos_edges[r_mat > 0 & p_mat < thresh] <- 1
  neg_edges[r_mat < 0 & p_mat < thresh] <- 1

  positive_edges <- which(pos_edges==1)
  negative_edges <- which(neg_edges==1)

  # get sum of all edges in TRAIN subs
  # divide by 2 to control since matrices are symmetric

  train_sumpos <- rep(0, N)
  train_sumneg <- rep(0, N)

  for (i in 1:N){
    train_sumpos[i] <- sum(train_array[,,i]*pos_edges, na.rm = TRUE)/2
    train_sumneg[i] <- sum(train_array[,,i]*neg_edges, na.rm = TRUE)/2
  }

  if(is.null(train_x) == TRUE){
    train_pos <- data.frame(sumpos=train_sumpos)
    train_neg <- data.frame(sumneg=train_sumneg)
    train_combined <- data.frame(sumpos=train_sumpos, sumneg=train_sumneg)
  }else{
    train_pos <- data.frame(sumpos=train_sumpos, x=train_x)
    train_neg <- data.frame(sumneg=train_sumneg, x=train_x)
    train_combined <- data.frame(sumpos=train_sumpos, sumneg=train_sumneg, x=train_x)
  }

  # build model on TRAIN subs

  if(edge=="separate"){

    if(fit == "linear"){
      if(sum(train_sumpos) != 0){
        fit_pos <- lm(train_behav ~ ., data = train_pos)
      }
      else{
        fit_pos <- NA
      }

      if(sum(train_sumneg) != 0){
        fit_neg <- lm(train_behav ~ ., data = train_neg)
      }
      else{
        fit_neg <- NA
      }
    }

    if(fit == "logistic"){
      if(sum(train_sumpos) != 0){
        fit_pos <- glm(train_behav ~ ., data = train_pos, family = "binomial")
      }

      else{
        fit_pos <- NA
      }

      if(sum(train_sumneg) != 0){
        fit_neg <- glm(train_behav ~ ., data = train_neg,  family = "binomial")
      }
      else{
        fit_neg <- NA
      }
    }

    if(fit == "multinom"){
      if(sum(train_sumpos) != 0){
        fit_pos <- multinom(train_behav ~ ., data = train_pos)
      }
      else{
        fit_pos <- NA
      }

      if(sum(train_sumneg) != 0){
        fit_neg <- multinom(train_behav ~ ., data = train_neg)
      }
      else{
        fit_neg <- NA
      }
    }

    return(list(r_mat=r_mat, p_mat=p_mat, positive_edges_matrix=pos_edges, negative_edges_matrix=neg_edges, positive_edges=positive_edges, negative_edges=negative_edges, positive_model=fit_pos, negative_model=fit_neg))

  }else if(edge=="combined"){

    if(fit == "linear"){
      if((sum(train_sumpos) != 0) | (sum(train_sumneg) != 0)){

        fit_combined <- lm(train_behav ~ ., data = train_combined)
      }
      else{
        fit_combined <- NA
      }
    }

    if(fit == "logistic"){
      if((sum(train_sumpos) != 0) | (sum(train_sumneg) != 0)){
        fit_combined <- glm(train_behav ~ ., data = train_combined, family = "binomial")
      }

      else{
        fit_combined <- NA
      }
    }

    if(fit == "multinom"){
      if((sum(train_sumpos) != 0) | (sum(train_sumneg) != 0)){

        fit_combined <- multinom(train_behav ~ ., data = train_combined)
      }
      else{
        fit_combined <- NA
      }
    }

    return(list(r_mat=r_mat, p_mat=p_mat, positive_edges_matrix=pos_edges, negative_edges_matrix=neg_edges, positive_edges=positive_edges, negative_edges=negative_edges, combined_model=fit_combined))

  }else{

    stop("model can only be fitted using either separate or combined edges")
  }



}


