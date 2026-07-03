#' collapseRDS
#'
#' @param rds.list list of rds objects
#'
#' @returns combined rds object
#' @export
collapseRDS <- function(rds.list) {
  comb <- rds.list[[1]]
  if (length(rds.list) > 1) {
    for (i in 2:length(rds.list)) {
      add <- rds.list[[i]]
      comb <- Biobase::combine(comb, add)
    }
    rds.obj <- comb
  } else {
    rds.obj <- comb
  }
  return(rds.obj)
}
