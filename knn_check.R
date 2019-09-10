setwd("~/Documents/MattesonResearch/AppliedBioPhysics/Paper2017")
load("datasets/gel_data.Rdata")
source("code/revisions-2/ECISfunctions.R")
library(e1071) # for "tune.knn"
library(class) # for "knn"

# Get data
freq.data <- gel.data

# Construct the frequency features
freqfeat <- freqFeat(freq.data)
r2h <- freqfeat[[1]]
mr <- freqfeat[[2]]
eof <- freqfeat[[3]]

# Extract key information from data - useful in analysis (see below)
n_freqs <- length(unique(freq.data$freq)) # number of unique frequencies
n_ts <- length(unique(freq.data$time_index)) # length of time series (associated with each observation and frequency)
n_ft <- n_freqs * n_ts # product of number of frequencies and length of time series
inds <- seq(1, nrow(freq.data), by = n_ft) # need one index associated with each well since Rb, alpha and cell type same across times and frequencies
rb <- freq.data$rb[inds]
a <- freq.data$alpha[inds]
cells <- freq.data$cell[inds]
freqs <- unique(freq.data$freq)
n_trials <- 20 

# Check KNN-CV
knn_check <- function(r2h, mr, eof, rb, a, freqs, cells, n_trials, n_feats){
  # Construct single/pairs/trios of features
  feature.data <- data.frame(r2h, mr, eof, rb, a) # r2h, mr, and eof matrices with columns corresponding to different frequencies
  meas_freq <- c(paste(rep(c("res2h", "maxres","endofrun"), each = length(freqs)), rep(freqs, 3)), "Rb", "Alpha")
  
  if (n_feats == 1){
    mf_length <- length(meas_freq)
    all_feats <- meas_freq
  } else if (n_feats == 2){
    feat_inds <- t(combn(1:length(meas_freq),2))
    mf_length <- nrow(feat_inds)
    all_feats <- paste(meas_freq[feat_inds[, 1]], meas_freq[feat_inds[, 2]])
  } else if (n_feats == 3){
    feat_inds <- t(combn(1:length(meas_freq),3))
    mf_length <- nrow(feat_inds)
    all_feats <- paste(meas_freq[feat_inds[, 1]], meas_freq[feat_inds[, 2]], meas_freq[feat_inds[, 3]])
  }
  class_rates <- data.frame(Features = all_feats, KNN = rep(0, mf_length), SE_KNN = rep(0, mf_length))
  cat("There are", mf_length, "features to go through.\n")
  
  for (k in 1:mf_length){
    # Construct data frame to be used in classification algorithms
    if (k %% 100 == 0){
      cat("k:", k, "\n")
    }
    if (n_feats == 1){
      feature1 <- feature.data[, k]
      f.data <- data.frame(cells, feature1)
    } else if (n_feats == 2){
      feature1 <- feature.data[, feat_inds[k, 1]]
      feature2 <- feature.data[, feat_inds[k, 2]]
      f.data <- data.frame(cells, feature1, feature2)
    } else if (n_feats == 3){
      feature1 <- feature.data[, feat_inds[k, 1]]
      feature2 <- feature.data[, feat_inds[k, 2]]
      feature3 <- feature.data[, feat_inds[k, 3]]
      f.data <- data.frame(cells, feature1, feature2, feature3)
    }
    
    # Create empty vectors to store out-of-sample predictive accuracy
    # for each of n_trials random splits of data using current feature and...
    kn <- rep(0, n_trials) # KNN-CV algorithm
    
    for (w in 1:n_trials){ # n_trials random splits of data
      set.seed(w) 
      
      # Get indicies of training set
      training_inds <- vector(length = length(cells))
      for (type in unique(cells)){
        inds <- which(cells == type)
        training_inds[sample(inds, 0.5*(length(inds)))] <- TRUE # randomly sample half of wells from each cell type
      }
      
      # Apply the three classification algorithms and out-of-sample predicitive accuracy
      knn.cross <- tune.knn(f.data[training_inds, -1], as.factor(f.data[training_inds, 1]), k = 1:20, tunecontrol=tune.control(sampling = "cross"), cross=10)
      knn.pred <- knn(as.matrix(f.data[training_inds, -1]), as.matrix(f.data[!training_inds, -1]), f.data[training_inds, 1], k=knn.cross$best.parameters[[1]]) # knn only works with matrices
      kn[w] <- mean(knn.pred == f.data[!training_inds, 1])
    }
    class_rates[k, 2] <- mean(kn)
    class_rates[k, 3] <- sd(kn)/sqrt(n_trials)
  }
  
  # KNN-CV in Column 2
  maxkn <- class_rates[which.max(class_rates[, 2]), ]
  
  return(maxkn)
}

# Apply function 
single_knn <- knn_check(r2h, mr, eof, rb, a, freqs, cells, n_trials, n_feats = 1)
single_knn # KNN-CV: 0.726 (compared to Tree: 0.714, LDA: 0.763, QDA: 0.762)
pair_knn <- knn_check(r2h, mr, eof, rb, a, freqs, cells, n_trials, n_feats = 2)
pair_knn # KNN-CV: 0.934 (compared to Tree: 0.930, LDA: 0.944, QDA: 0.951, RDA: 0.969)
trio_knn <- knn_check(r2h, mr, eof, rb, a, freqs, cells, n_trials, n_feats = 3)
trio_knn # KNN-CV: 0.982 (compared to Tree: 0.970, LDA: 0.981, QDA: 0.978, RDA: 0.992)

knn_check <- list(single_knn = single_knn, pair_knn = pair_knn, trio_knn = trio_knn)

save(knn_check, file = "/Users/megangelsinger/Documents/MattesonResearch/AppliedBioPhysics/Paper2017/InternationalJournalOfBiostatistics/Round_2_Minor_Revisions/figs_tables/knn_check.Rdata")
 