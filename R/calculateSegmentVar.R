calculateSegmentVar <- function(abs_seg = NULL, abs_cn = NULL) {
  segRLE <- rle(as.numeric(abs_seg))

  segVar <- c()
  for (i in 1:length(segRLE$lengths)) {
    if (i == 1) {
      strt_idx <- 1
      end_idx <- segRLE$lengths[i]
      segVar <- append(segVar, var(abs_cn[strt_idx:end_idx]))
    } else {
      start_idx <- max(cumsum(segRLE$lengths[1:i - 1]))
      strt_idx <- start_idx + 1
      if (i == length(segRLE$lengths)) {
        end_idx <- segRLE$lengths[i] + strt_idx - 1
      } else {
        end_idx <- segRLE$lengths[i] + strt_idx
      }
      segVar <- append(segVar, var(abs_cn[strt_idx:end_idx]))
    }
  }

  medianVar <- median(segVar)
  return(medianVar)
}
