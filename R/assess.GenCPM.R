#' Given a summary performance measures for the GenCPM model
#'
#' @import psych
#' @import pROC
#' @param object Returned GenCPM object from `.GenCPM` or `.regularized.GenCPM` functions. 
#' @param model A character string representing one of the built-in regression models. 
#' “linear” for `linear.GenCPM` and `linear.regularized.GenCPM`; 
#' “logistic” for `logit.GenCPM` and `logit.regularized.GenCPM`;
#' “multinom” for `multinom.GenCPM` and `multinom.regularized.GenCPM`; 
#' and  “cox” for `cox.GenCPM` and `cox.regularized.GenCPM`. 
#' The default is “linear”.
#' @param edge Usage of edges to fit models, and it should be decided by the edge usage in the "object" input.
#' “seperate” for fitting two separate models using positive edges 
#' and negative edges respectively, and “combined” for fitting only one model 
#' use all edges selected. The default is “separate”.

#' @return A list contains metrics assessing the model performance (MSE, AUC, C-index, etc), predicted response and actual response.
#' @export

assess.GenCPM <- function(object, model="linear", edge="separate"){
  
  if(!(model %in% c("linear", "logistic", "cox", "multinom"))){
    stop("model must be one of the specific types")
  }
  
  if(!(edge %in% c("separate", "combined"))){
    stop("edge must be separate or combined")
  }
  

  if(model=="linear"){
    
    if(edge == "separate"){
      
      behav_actual <- as.vector(object$actual_behavior)
      N <- length(behav_actual)
      behav_pred_pos <- as.vector(object$positive_predicted_behavior)
      behav_pred_neg <- as.vector(object$negative_predicted_behavior)
      
      MSE_pos <- sum((behav_pred_pos-behav_actual)^2, na.rm = TRUE)/N
      corr_pos <- corr.test(x = as.vector(behav_pred_pos), y = as.vector(behav_actual), adjust = "none", ci=F)
      R_pos <- corr_pos$r

      MSE_neg <- sum((behav_pred_neg-behav_actual)^2, na.rm = TRUE)/N
      corr_neg <- corr.test(x = as.vector(behav_pred_neg), y = as.vector(behav_actual), adjust = "none", ci=F)
      R_neg <- corr_neg$r

      return(list(positive_predicted_behavior=behav_pred_pos, 
                  negative_predicted_behavior=behav_pred_neg, 
                  actual_behavior=behav_actual, 
                  positive_r=R_pos, negative_r=R_neg, 
                  positive_MSE=MSE_pos, negative_MSE=MSE_neg))
      
    }else if(edge=="combined"){
      
      behav_actual <- as.vector(object$actual_behavior)
      N <- length(behav_actual)
      behav_pred <- as.vector(object$predicted_behavior)

      MSE <- sum((behav_pred-behav_actual)^2, na.rm = TRUE)/N
      corr <- corr.test(x = as.vector(behav_pred), y = as.vector(behav_actual), adjust = "none", ci=F)
      R <- corr$r

      return(list(predicted_behavior=behav_pred, 
                  actual_behavior=behav_actual, 
                  r=R, MSE=MSE))
      
    }
    
  }
  
  if(model=="logistic"){
      
      
      if(edge == "separate"){
        
        behav_actual <- as.vector(object$actual_behavior)
        N <- length(behav_actual)
        behav_pred_pos <- as.vector(object$positive_predicted_behavior)
        behav_pred_neg <- as.vector(object$negative_predicted_behavior)
        
        AUC_pos <- auc(behav_actual, behav_pred_pos)
        AUC_neg <- auc(behav_actual, behav_pred_neg)
        
        return(list(positive_predicted_behavior=behav_pred_pos, 
                    negative_predicted_behavior=behav_pred_neg, 
                    actual_behavior=behav_actual, 
                    positive_AUC=AUC_pos, negative_AUC=AUC_neg))
        
      }else if(edge == "combined"){
        
        behav_actual <- as.vector(object$actual_behavior)
        N <- length(behav_actual)
        behav_pred <- as.vector(object$predicted_behavior)
        
        AUC <- auc(behav_actual, behav_pred)
        
        return(list(predicted_behavior=behav_pred, 
                    actual_behavior=behav_actual, 
                    AUC=AUC))
        
      }
      
    
  }
  
  if(model=="multinom"){
    
    if(edge == "separate"){

      behav_actual <- as.vector(object$actual_behavior)
      N <- length(behav_actual)
      behav_pred_pos <- as.numeric(as.vector(object$positive_predicted_behavior))
      behav_pred_neg <- as.numeric(as.vector(object$negative_predicted_behavior))
      
      AUC_pos <- multiclass.roc(behav_actual, behav_pred_pos)$auc
      AUC_neg <- multiclass.roc(behav_actual, behav_pred_neg)$auc
      
      return(list(positive_predicted_behavior=behav_pred_pos, 
                  negative_predicted_behavior=behav_pred_neg, 
                  actual_behavior=behav_actual, 
                  positive_AUC=AUC_pos, negative_AUC=AUC_neg))
      
    }else if(edge == "combined"){
      
      behav_actual <- as.vector(object$actual_behavior)
      N <- length(behav_actual)
      behav_pred <- as.numeric(as.vector(object$predicted_behavior))
      
      AUC <- multiclass.roc(behav_actual, behav_pred)$auc
      
      return(list(predicted_behavior=behav_pred, 
                  actual_behavior=behav_actual, 
                  AUC=AUC))
      
    }
    
  }
  
  if(model=="cox"){
    
    if(edge=="separate"){
      
      pos_lp <- as.vector(object$positive_predicted_linear_predictor)
      neg_lp <- as.vector(object$negative_predicted_linear_predictor)
      
      status <- as.vector(object$actual_status)
      time <- as.vector(object$actual_time)

      ts <-  data.frame("time" = time, "status" = status)
      pos_Cindex <- Cindex(pos_lp, ts, weights = rep(1, nrow(ts)))
      neg_Cindex <- Cindex(neg_lp, ts, weights = rep(1, nrow(ts)))

      return(list(positive_predicted_linear_predictor=pos_lp, 
                  negative_predicted_linear_predictor=neg_lp,
                  actual_status=status,
                  actual_time=time, 
                  positive_cindex=pos_Cindex,
                  negative_cindex=neg_Cindex))
      
    }else if(edge=="combined"){
      
      lp <- as.vector(object$predicted_linear_predictor)
      status <- as.vector(object$actual_status)
      time <- as.vector(object$actual_time)
      
      ts <-  data.frame("time" = time, "status" = status)
      Cindex <- Cindex(lp, ts, weights = rep(1, nrow(ts)))

      return(list(predicted_linear_predictor=lp, 
                  actual_status=status,
                  actual_time=time, cindex=Cindex))
    }
    
  }
  
}