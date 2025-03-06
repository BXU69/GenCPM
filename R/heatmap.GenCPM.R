#' Generate heatmaps for the GenCPM object based on 10-node network label in Shen268 atlas.
#' @import ComplexHeatmap
#' @import circlize
#' @param cpm Returned GenCPM object from `.GenCPM` or `.regularized.GenCPM` functions. 
#' @param foldThreshold the edges selected for over this many folds will be plotted. If set to .5, the edges selected at least half of the time are plotted.  
#' @export

 
heatmap.GenCPM <- function(cpm, foldThreshold = .5){
  
  nfold <- length(cpm$positive_edges)
  
  pos_number <- unlist(cpm$positive_edges)
  neg_number <- unlist(cpm$negative_edges)
  pos_count <- as.data.frame(table(pos_number))
  neg_count <- as.data.frame(table(neg_number))
  
  # add the threshold criteria here
  pos_count <- pos_count[pos_count$Freq >= foldThreshold*nfold,]
  neg_count <- neg_count[neg_count$Freq >= foldThreshold*nfold,]
  
  if(nrow(pos_count) == 0 && nrow(neg_count) == 0){
    
    return("there is no edge selected at this foldThresh.")
    
  }else{
    
    pos_prop <- pos_count$Freq/sum(pos_count$Freq)
    neg_prop <- neg_count$Freq/sum(neg_count$Freq)
    
    pos_count$pos_number <- as.numeric(as.character(pos_count$pos_number))
    neg_count$neg_number <- as.numeric(as.character(neg_count$neg_number))
    
    # calculate the postion of each edge
    
    pos_position1 = ceiling(pos_count$pos_number / 268) 
    pos_position2 = pos_count$pos_number - (ceiling(pos_count$pos_number / 268) - 1) * 268
    
    neg_position1 = ceiling(neg_count$neg_number / 268) 
    neg_position2 = neg_count$neg_number - (ceiling(neg_count$neg_number / 268) - 1) * 268
    
    # load the shen_268_network_labels and 10-node network labels
    
    # load("shen_268_network_labels.RData")
    node268 <- shen_268_network_labels$`10-node_network`
    
    # match position with 10-node labels
    
    pos_position <- data.frame("position1"=node268[pos_position1], "position2"=node268[pos_position2], "proportion"=pos_prop)
    neg_position <- data.frame("position1"=node268[neg_position1], "position2"=node268[neg_position2], "proportion"=neg_prop)
    
    # set the label of 10-node network
    
    # label <- c("Medial Frontal", "Fronto-parietal", "Default Mode", "Motor", 
    # "Visual I", "Visual II", "Visual Association", "Limbic", "Basal Ganglia", "Cerebellum")
    label <- c("MG","FP","DM","MOT","V1","V2","VA","LIM","BG","CER")
    
    pos_position$position1 <- label[pos_position$position1]
    pos_position$position2 <- label[pos_position$position2]
    
    neg_position$position1 <- label[neg_position$position1]
    neg_position$position2 <- label[neg_position$position2]
    
    # get the upper triangle part of the frequancy matrix
    # for positive edges
    
    pp <- data.frame(matrix(ncol = 10, nrow = 10))
    colnames(pp) <- label
    row.names(pp) <- label
    
    for (i in 1:10){
      for (j in 1:10){
        ppe <- sum(pos_position$proportion[which(pos_position$position1==label[i]&pos_position$position2==label[j])])
        if(length(ppe)!=0){
          pp[i,j] <- ppe
        }else{
          pp[i,j] <- 0
        }
      }
    }
    
    pp <- data.matrix(t(pp)+pp)
    diag(pp) <- diag(pp)/2
    
    # for negative edges
    np <- data.frame(matrix(ncol = 10, nrow = 10))
    colnames(np) <- label
    row.names(np) <- label
    
    for (i in 1:10){
      for (j in 1:10){
        npe <- sum(neg_position$proportion[which(neg_position$position1==label[i]&neg_position$position2==label[j])])
        if(length(npe)!=0){
          np[i,j] <- npe
        }else{
          np[i,j] <- 0
        }
      }
    }
    
    np <- data.matrix(np+t(np))
    diag(np) <- diag(np)/2
    np <- -np
    
    # draw the heatmap
    max_pos <- max(unique(pp))
    max_neg <- abs(min(unique(np)))
    max <- max(max_pos, max_neg)
    col1 = colorRamp2(c(-max, 0, max), c("blue","white", "red"))
    
    ht1 <- Heatmap(pp, rect_gp = gpar(type = "none"), col=col1,
                   width = 10*unit(10, "mm"),
                   height = 10*unit(10, "mm"),
                   heatmap_legend_param = list(title=" ",at = c(-max, max), 
                                               labels = c("negative", "positive"),
                                               direction = "horizontal",
                                               legend_width = unit(8, "cm")),
                   cluster_rows = FALSE, cluster_columns = FALSE, 
                   row_names_side = "left",
                   column_names_rot = 45,
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(i >= j) {
                       grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "dimgray"))
                     }
                   })
    
    ht2 <- Heatmap(np, rect_gp = gpar(type = "none"), col=col1,  
                   width = 10*unit(10, "mm"),
                   height = 10*unit(10, "mm"),
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   column_names_side = "top",
                   show_heatmap_legend = FALSE,
                   column_names_rot = 45,
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(i <= j) {
                       grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "dimgray"))
                     }
                   }
    )
    
    
    draw(ht1 + ht2, ht_gap = unit(-70, "mm"), heatmap_legend_side = "bot")
    
    
  }
  
}
