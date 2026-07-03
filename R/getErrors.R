getErrors <- function(data=NULL){
  # clonality is a legacy name for MAE
  clonality <- mean(abs(data$errors))
  # Root Mean Squared Error
  rmse <- sqrt(mean(data$errors^2))
  MedianSegVar <- calculateSegmentVar(abs_seg = data$abs_seg,
                                      abs_cn = data$abs_cn)
  return(c(clonality=clonality,
           rmse=rmse,
           MedianSegVar=MedianSegVar))
}
