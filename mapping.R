################# Approach 1 #################
# Load required libraries
library(vcfR)
library(ggplot2)
library(qtl)

cross <- read.cross("csv", "", "Approach1.csv", 
                    estimate.map = TRUE, map.function = "kosambi")

# Convert to RIL population and jitter overlapping markers
cross <- convert2riself(cross)
cross <- jittermap(cross)

# Estimate initial map
map <- est.map(cross, map.function = "kosambi", error.prob = 0.001)
plotMap(map)
summary.map(map)


marker_info <- function(map_in, dist_in) {
  temp_out <- numeric()
  gt <- geno.table(map_in, scanone.output = TRUE)
  for (i in unique(gt$chr)) {
    temp_chr <- gt[gt$chr == i, ]
    chr_out <- calc_distance(temp_chr, dist_in)
    if (is.data.frame(chr_out)) {
      temp_out <- rbind(temp_out, chr_out)
    }
  }
  return(temp_out)
}

calc_distance <- function(chr_in, dist_in) {
  temp_out <- numeric()
  for (i in 2:(nrow(chr_in) - 1)) {
    pre_dist <- chr_in[i, 'pos'] - chr_in[i - 1, 'pos']
    post_dist <- chr_in[i + 1, 'pos'] - chr_in[i, 'pos']
    if (pre_dist > dist_in && post_dist > dist_in) {
      temp_out <- rbind(temp_out, chr_in[i, ])
    }
  }
  return(temp_out)
}

thinning_loop <- function(map_in) {
  gt <- geno.table(map_in, scanone.output = TRUE)
  snps <- rownames(gt)
  out <- data.frame(id = character(), dist = numeric(), pos = numeric(), 
                    neglog10P = numeric(), missing = numeric())
  for (i in snps[2:length(snps)]) {
    t <- which(rownames(gt) == i)
    pre <- t - 1
    if (round(gt$pos[t]) == round(gt$pos[pre])) {
      dist <- round(gt$pos[t]) - round(gt$pos[pre])
      best <- if (gt$neglog10P[t] > gt$neglog10P[pre]) gt[t, ] else gt[pre, ]
      best$id <- rownames(best)
      out <- rbind(out, data.frame(id = best$id, dist = dist, pos = best$pos, 
                                   neglog10P = best$neglog10P, missing = best$missing))
    }
  }
  return(unique(out$id))
}

calc_distance2 <- function(chr_in, dist_in, num_mar_in, inner_dist) {
  temp_out <- numeric()
  window_width <- num_mar_in - 1
  for (i in 2:(nrow(chr_in) - num_mar_in)) {
    if (i + num_mar_in < nrow(chr_in)) {
      window_range <- chr_in[(i + window_width), 'pos'] - chr_in[i, 'pos']
      if (window_range < inner_dist) {
        pre_dist <- chr_in[i, 'pos'] - chr_in[i - 1, 'pos']
        post_dist <- chr_in[i + num_mar_in, 'pos'] - chr_in[i + window_width, 'pos']
        if (pre_dist > dist_in && post_dist > dist_in) {
          temp_out <- rbind(temp_out, chr_in[i:(i + window_width), ])
        }
      }
    }
  }
  return(temp_out)
}

marker_info2 <- function(map_in, dist_in, num_mar_in, inner_dist) {
  temp_out <- numeric()
  gt <- geno.table(map_in, scanone.output = TRUE)
  for (i in unique(gt$chr)) {
    temp_chr <- gt[gt$chr == i, ]
    if (nrow(temp_chr) > num_mar_in) {
      chr_out <- calc_distance2(temp_chr, dist_in, num_mar_in, inner_dist)
      if (is.data.frame(chr_out)) {
        temp_out <- rbind(temp_out, chr_out)
      }
    }
  }
  if (anyNA(temp_out)) {
    temp_out <- na.omit(temp_out)
  }
  return(temp_out)
}

distbetweenMarkers <- function(data_in) {
  data_in$dist <- NA
  colnames(data_in)[1] <- 'pos'
  for (i in 2:nrow(data_in)) {
    data_in[i, 'dist'] <- round(data_in[i, 'pos'] - data_in[i - 1, 'pos'], 2)
  }
  return(data_in)
}

bmarkers <- thinning_loop(cross)

# Remove co-localized markers and re-estimate the map
cross2 <- drop.markers(cross, bmarkers)
nm <- est.map(cross2, map.function = "kosambi", error.prob = 0.001)
cross2 <- replace.map(cross2, nm)


# Plot recombination frequencies
plotRF(cross2, mark.diagonal = TRUE, col.scheme = 'redblue')
plotRF(cross2, chr = '1', mark.diagonal = TRUE, col.scheme = 'redblue')
abline(v = seq(0, 900, 10), lwd = 1.7)



################# Approach 2 #################
# Load in the dataset (formatted for R/qtl, created via vcf_to_rqtl.R)
cross <- read.cross("csv", "", "Approach2.csv", 
                    estimate.map = TRUE, map.function = "kosambi")

# Convert to recombinant inbred lines
cross <- convert2riself(cross)

# Data summary and initial processing
summary(cross)
cross <- jittermap(cross)

# Estimate recombination frequencies and plot initial map
map <- est.map(cross, map.function = "kosambi", error.prob = 0.001)
plotMap(map)

# Example: Recombination frequency plot for chromosome 17 to identify erroneous markers
plotRF(cross, chr = "17", mark.diagonal = TRUE, col.scheme = "redblue")
abline(v = seq(0, 900, 10), lwd = 1.7)

# Drop problematic or redundant markers 
drop <- c(
  "c3m176307608", "c6m11862214", "c9m127843085",
  "c10m49576888", "c10m49576917", "c10m142996839",
  "c10m190339292", "c10m190339301", "c10m190339331",
  "c10m190339343", "c10m190339358", "c10m190469822",
  "c10m190469837", "c13m24509319", "c13m175660374",
  "c17m195518598", "c17m195518737"
)

cross2 <- drop.markers(cross, drop)

# Re-estimate genetic map and visualize
nm <- est.map(cross2, map.function = "kosambi", error.prob = 0.001)
plotMap(nm)
plotRF(cross2, mark.diagonal = TRUE, col.scheme = "redblue")
summary(cross2)

# Heatmap of LOD and recombination values
heatMap(cross2, what = c("both", "lod", "rf"), lmax = 12, rmin = 0,
        markDiagonal = TRUE, color = rev(rainbow(256, start = 0, end = 2/3)))

# Calculate genotype probabilities
cross2 <- calc.genoprob(cross2, step = 1, map.function = "kosambi", error.prob = 0.001)
summaryMap(cross2)
# 1D
oned <- scanone(cross2, pheno.col = 1, model = "normal", method = "hk")
plot(oned, main = "Interval Mapping Results")
summary(oned)
#CIM
cim.out <- cim(cross2, pheno.col = 1, map.function = "kosambi")
plot(cim.out, col = "darkblue", main = "CIM-QTL Results for Nectar Volume", ylab = "LOD score")
abline(h = 1.5, col = "grey", lty = 2)
summary(cim.out)
max(cim.out)

# Permutation test for CIM
cim.perm <- cim(cross2, pheno.col = 1, n.perm = 1000, map.function = "kosambi")
summary(cim.perm, alpha = 0.05)

out2d <- scantwo(cross2, pheno.col = 1, model = "normal", method = "hk")
plot(out2d, point.at.max = TRUE, main = "2D Genome Scan LOD Scores")
summary(out2d)
max(out2d)

# Permutation test for 2D scan
operm2 <- scantwo(cross2, method = "hk", n.perm = 1000)
summary(operm2, alpha = 0.05)

################# Approach 3 #################
### 1. Load and Prepare Data

cross <- read.cross("csv", "", "Approach3.csv", 
                    estimate.map = TRUE, map.function = "kosambi")

summary(cross)
cross <- jittermap(cross)
map <- est.map(cross, map.function = "kosambi", error.prob = 0.001)
cross <- replace.map(cross, map)

### 2. Marker Thinning Functions
thinning_loop <- function(map_in) {
  gt <- geno.table(map_in, scanone.output = TRUE)
  snps <- rownames(gt)
  out <- data.frame()
  for (i in 2:length(snps)) {
    pre <- i - 1
    if (round(gt$pos[i]) == round(gt$pos[pre])) {
      selected <- if (gt$neglog10P[i] > gt$neglog10P[pre]) gt[i, ] else gt[pre, ]
      out <- rbind(out, data.frame(id = rownames(selected), 
                                   dist = gt$pos[i] - gt$pos[pre],
                                   pos = selected$pos,
                                   neglog10P = selected$neglog10P,
                                   missing = selected$missing))
    }
  }
  return(unique(out$id))
}


### 3. Drop Erroneous SNPs Based on RF Maps and Frequencies
# --- Round 1 thinning ---
badmarkers <- thinning_loop(cross)
cross2 <- drop.markers(cross, badmarkers)
nm <- est.map(cross2, error.prob = 0.001, map.function = "kosambi")
cross2 <- replace.map(cross2, nm)

# --- Round 2 drop ---
drop <- c("c2m25220136","c2m29398816","c2m67080281","c2m139916392","c2m177964508",
          "c3m43133858","c3m109484596","c3m141613925","c3m171291089",
          "c5m20952397","c5m26395632","c5m36156808","c5m37628308","c5m40297745",
          "c5m91887832","c5m92720616","c5m112048741","c5m113611894","c5m115373892",
          "c5m119174365","c5m121238655","c5m121470413","c5m125751016","c5m127304124",
          "c5m127709570","c5m133714784","c5m147020312","c5m145130924","c5m148052011",
          "c5m149588749","c5m161843721","c5m167939353","c5m169414598","c5m173800973",
          "c6m1241403","c6m12165701","c7m63896184","c7m97584345","c7m131408436",
          "c11m10028095","c11m49925322","c11m76160435","c11m101864327","c11m116417706",
          "c11m155343142","c11m157438671","c11m162343611","c11m176967372",
          "c12m4456796","c12m16091016","c12m161567797","c13m10667921","c13m17096571",
          "c13m41677980","c14m178101413","c14m178383270","c14m178988886",
          "c14m180093446","c14m182683330","c15m4596433","c15m12698266","c15m15700937",
          "c15m17279131","c15m26678562","c15m53717050","c15m60340806","c15m111257309",
          "c15m133720182","c15m156915216","c15m157279313","c15m159428618",
          "c16m1944189","c16m90482469","c16m96393690","c16m165248861",
          "c16m185325020","c16m215426458","c17m12366982")
cross3 <- drop.markers(cross2, drop)
nm2 <- est.map(cross3, error.prob = 0.001, map.function = "kosambi")
cross3 <- replace.map(cross3, nm2)

# --- Round 3 drop ---
drop2 <- c("c5m37954255","c8m62200223","c9m129499268","c11m179907155","c15m149713068",
           "c15m149363120", "c5m20078994", "c5m21314710", "c5m24310267", "c5m25934797",
           "c5m26847421", "c5m27655961", "c5m34378687", "c5m35733592", "c5m93400760",
           "c5m94509439", "c5m95401964", "c5m95856887", "c5m100265586", "c5m105958548",
           "c5m109914270", "c5m112760010", "c5m114271535", "c5m120061728", "c5m125465003",
           "c5m130508843", "c5m140661729", "c5m146567121", "c5m147581970", "c5m148681639",
           "c5m148871649", "c5m152119174", "c5m155797506", "c5m157223236", "c5m158139733",
           "c5m158635081", "c5m159359941", "c5m160182037", "c5m162399568", "c5m169101560")
cross4 <- drop.markers(cross3, drop2)
nm3 <- est.map(cross4, error.prob = 0.001, map.function = "kosambi")
cross4 <- replace.map(cross4, nm3)

### 4. QTL Mapping & Model Fitting
# Calculate genotype probabilities
cross4 <- calc.genoprob(cross4, step = 1, map.function = "kosambi", error.prob = 0.001)

# Interval Mapping
oned <- scanone(cross4, pheno.col = 1, model = "normal", method = "hk")
plot(oned, main = "Interval Mapping Results")
summary(oned)

# Composite Interval Mapping
test_cim2 <- cim(cross4, pheno.col = 1, map.function = "kosambi", method = "hk",
                 window = 100, n.marcovar = 2)
summary(test_cim2)
plot(test_cim2, col = "darkblue", main = "CIM-QTL Results for Nectar Volume", ylab = "LOD score")
abline(h = 3, col = "grey", lty = 2)
write.csv(test_cim2, "cim_results.csv", row.names = FALSE)

# Refine QTL model
qtl.model2 <- makeqtl(cross4,
                      chr = c(10, 11, 15, 16, "1", "1", "6", "9", 10),
                      pos = c(402.9418, 3.211002, 16.215, 32.94041, 88.022, 191.022,
                              4.987652, 248.8825, 293.635),
                      what = "prob")
fit.result <- fitqtl(cross4,
                     qtl = qtl.model2,
                     formula = y ~ Q1:Q2 + Q3:Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q4,
                     method = "hk", dropone = TRUE, get.ests = TRUE)
summary(fit.result)

# LOD Interval Estimation
rqtl <- refineqtl(cross4,
                  pheno.col = 1,
                  model = "normal",
                  qtl = qtl.model2,
                  formula = y ~ Q1:Q2 + Q3:Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q1:Q4)
plotLodProfile(rqtl)

# LOD intervals
for (i in 1:9) {
  print(lodint(rqtl, qtl.index = i, drop = 1.5))
}

