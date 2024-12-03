# Author: Ashley Barstow

# Run thinning Loop to 
# Define the marker_info function
marker_info <- function(map_in, dist_in) {
  temp_out <- numeric()
  gt <- geno.table(map_in, scanone.output = T)
  for (i in unique(gt$chr)) {
    temp_chr <- gt[gt$chr == i, ]
    chr_out <- calc_distance(temp_chr, dist_in)
    if (is.data.frame(chr_out) == TRUE) {
      temp_out <- rbind.data.frame(temp_out, chr_out)
    }
  }
  return(temp_out)
}

# Define the calc_distance function
calc_distance <- function(chr_in, dist_in) {
  temp_out <- numeric()
  temp_length <- nrow(chr_in)
  for (i in seq(2, (temp_length - 1), 1)) {
    temp_dist_pre <- chr_in[i, 'pos'] - chr_in[i - 1, 'pos']
    temp_dist_post <- chr_in[i + 1, 'pos'] - chr_in[i, 'pos']
    if (temp_dist_post > dist_in & temp_dist_pre > dist_in) {
      temp_out2 <- chr_in[i,]
      temp_out <- rbind.data.frame(temp_out, temp_out2)
    }
  }
  return(temp_out)
}

# Define the thinning_loop function
thinning_loop <- function(map_in) {
  gt <- geno.table(map_in, scanone.output = TRUE)
  snp <- as.character(rownames(gt))
  out <- data.frame(id = as.character(), dist = as.numeric(), pos = as.numeric(), neglog10P = as.numeric(), missing = as.numeric())
  len <- as.numeric(length(snp))
  for (i in snp[2:len]) {
    t <- which(rownames(gt) == i)
    pre <- t - 1
    if (round(gt$pos[t], digits = 0) == round(gt$pos[pre], digits = 0)) {
      dist <- round(gt$pos[t], digits = 0) - round(gt$pos[pre], digits = 0)
      if (gt$neglog10P[t] > gt$neglog10P[pre]) {
        ugh <- gt[t, c(2:4)]
      } else {
        ugh <- gt[pre, c(2:4)]
      }
      ugh[, "id"] <- paste(rownames(ugh))
      out <- rbind(out, cbind(id = ugh[[4]], dist = dist, pos = ugh[[1]], neglog10P = ugh[[2]], missing = ugh[[3]]))
    }
  }
  todrop <- unique(out$id)
  return(todrop)
}

# Define the calc_distance2 function
calc_distance2 <- function(chr_in, dist_in, num_mar_in, inner_dist) {
  temp_out <- numeric()
  temp_length <- nrow(chr_in)
  window_width <- num_mar_in - 1
  for (i in seq(2, (temp_length - num_mar_in), 1)) {
    if (i + num_mar_in < nrow(chr_in)) {
      if (chr_in[(i + window_width), 'pos'] - chr_in[i, 'pos'] < inner_dist) {
        temp_dist1 <- chr_in[i, 'pos'] - chr_in[i - 1, 'pos']
        temp_dist2 <- chr_in[i + num_mar_in, 'pos'] - chr_in[i + window_width, 'pos']
        if (temp_dist1 > dist_in & temp_dist2 > dist_in) {
          temp_out2 <- chr_in[i:(i + (num_mar_in - 1)), ]
          temp_out <- rbind.data.frame(temp_out, temp_out2)
        }
      }
    }
  }
  return(temp_out)
}

# Define the marker_info2 function
marker_info2 <- function(map_in, dist_in, num_mar_in, inner_dist) {
  temp_out <- numeric()
  gt <- geno.table(map_in, scanone.output = T)
  for (i in unique(gt$chr)) {
    temp_chr <- gt[gt$chr == i, ]
    if (nrow(temp_chr) > num_mar_in) {
      chr_out <- calc_distance2(temp_chr, dist_in, num_mar_in, inner_dist)
      if (is.data.frame(chr_out) == TRUE) {
        temp_out <- rbind.data.frame(temp_out, chr_out, inner_dist)
      }
    }
  }
  if (sum(is.na(temp_out)) != 0) {
    temp_out <- na.omit(temp_out)
  }
  return(temp_out)
}

# Define the distbetweenMarkers function
distbetweenMarkers <- function(data_in) {
  data_in$dist <- NA
  colnames(data_in)[1] <- 'pos'
  for (i in 2:nrow(data_in)) {
    data_in[i, 'dist'] <- round(data_in[i, 'pos'] - data_in[i - 1, 'pos'], 2)
  }
  return(data_in)
}


badmarkers <- thinning_loop(cross)
