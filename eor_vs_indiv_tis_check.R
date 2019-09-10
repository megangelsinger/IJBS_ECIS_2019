load("gel_data.Rdata")
source("ECISfunctions.R")
library(tree)  # for 'tree'
library(MASS) # for 'lda' and 'qda'

############################## Get data ##################################
freq.data <- gel.data

# Construct the frequency features (using last five time indices individually, as opposed to averaged through EOR)
freqfeat <- freqFeat(freq.data)
r2h <- freqfeat[[1]]
mr <- freqfeat[[2]]
ti60 <- matrix(freq.data[freq.data$time_index == 60, "res"], nrow = 210, ncol = 9, byrow = TRUE)
ti61 <- matrix(freq.data[freq.data$time_index == 61, "res"], nrow = 210, ncol = 9, byrow = TRUE)
ti62 <- matrix(freq.data[freq.data$time_index == 62, "res"], nrow = 210, ncol = 9, byrow = TRUE)
ti63 <- matrix(freq.data[freq.data$time_index == 63, "res"], nrow = 210, ncol = 9, byrow = TRUE)
ti64 <- matrix(freq.data[freq.data$time_index == 64, "res"], nrow = 210, ncol = 9, byrow = TRUE)

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

########################### Check with old features (minus EOR) and new features #################
check_tis <- function(r2h, mr, ti60, ti61, ti62, ti63, ti64, rb, a, freqs, cells, n_trials, n_feats){
  # Construct single/pairs/trios of features
  feature.data <- data.frame(r2h, mr, ti60, ti61, ti62, ti63, ti64, rb, a) # r2h, mr, and ti60-604 matrices with columns corresponding to different frequencies
  meas_freq <- c(paste(rep(c("res2h", "maxres","ti60", "ti61", "ti62", "ti63", "ti64"), each = length(freqs)), rep(freqs, 7)), "Rb", "Alpha")
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
  class_rates <- data.frame(Features = all_feats, Tree = rep(0, mf_length), SE_Tree = rep(0, mf_length), LDA = rep(0, mf_length), SE_LDA = rep(0, mf_length), QDA = rep(0, mf_length), SE_QDA = rep(0, mf_length))
  cat("There are", mf_length, "features to go through.\n")
  
  for (k in 1:mf_length){
    if (k %% 100 == 0){
      cat("k:", k, "\n")
    }
    # Construct data frame to be used in classification algorithms
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
    tr <- rep(0, n_trials) # classification tree algorithm
    l <- rep(0, n_trials) # LDA
    q <- rep(0, n_trials) # QDA
    
    for (w in 1:n_trials){ # n_trials random splits of data
      set.seed(w) 
      
      # Get indicies of training set
      training_inds <- vector(length = length(cells))
      for (type in unique(cells)){
        inds <- which(cells == type)
        training_inds[sample(inds, 0.5*(length(inds)))] <- TRUE # randomly sample half of wells from each cell type
      }
      
      # Apply the three classification algorithms and out-of-sample predicitive accuracy
      tree.fit <- tree(cells ~ ., data = f.data, subset = training_inds)
      lda.fit <- lda(cells ~ ., data = f.data, subset = training_inds)
      qda.fit <- qda(cells ~ ., data = f.data, subset = training_inds)
      
      tree.pred <- predict(tree.fit, f.data[!training_inds, ],type = "class")
      tree.class <- tree.pred
      tr[w] <- mean(tree.class == f.data[!training_inds, 1])
      
      lda.pred <- predict(lda.fit, f.data[!training_inds, ])
      lda.class <- lda.pred$class
      l[w] <- mean(lda.class == f.data[!training_inds, 1])
      
      qda.pred <- predict(qda.fit, f.data[!training_inds, ])
      qda.class <- qda.pred$class
      q[w] <- mean(qda.class == f.data[!training_inds, 1])
      
    }
    class_rates[k, 2] <- mean(tr)
    class_rates[k, 3] <- sd(tr)/sqrt(n_trials)
    class_rates[k, 4] <- mean(l)
    class_rates[k, 5] <- sd(l)/sqrt(n_trials)
    class_rates[k, 6] <- mean(q)
    class_rates[k, 7] <- sd(q)/sqrt(n_trials)
  }
  
  # Tree in Column 2
  maxt <- class_rates[which.max(class_rates[, 2]), ]
  # LDA in Column 4
  maxl <- class_rates[which.max(class_rates[, 4]), ]
  # QDA in Column 6
  maxq <- class_rates[which.max(class_rates[, 6]), ]
  # Combine
  best_feats <- rbind(maxt, maxl, maxq)
  
  # Return data associated with the best features
  best_feat_data <- array(dim = c(nrow(r2h), n_feats, 3)) # for each classification method, matrix of size n_wells x n_feats
  for (i in 1:3){
    names_best_feats <- unlist(strsplit(as.character(best_feats[i, 1]), split = " "))
    rb_ind <- which(names_best_feats == "Rb")
    alpha_ind <- which(names_best_feats == "Alpha")
    interum_mat <- matrix(NA, nrow = nrow(r2h), ncol = length(names_best_feats)) # must accommodate Rb and Alpha which don't have frequencies, still want features in same order
    if (length(rb_ind) != 0){
      interum_mat[, rb_ind] <- feature.data[[which(meas_freq == "Rb")]]
    } 
    if (length(alpha_ind) != 0){
      interum_mat[, alpha_ind] <- feature.data[[which(meas_freq == "Alpha")]]
    }
    other_inds <- which(names_best_feats != "Rb"& names_best_feats != "Alpha")
    for (j in other_inds[seq(1, length(other_inds), by = 2)]){
      interum_mat[, j] <- feature.data[[which(meas_freq == paste(names_best_feats[j], names_best_feats[j+1]))]]
    }
    best_feat_data[,,i] <- interum_mat[, colSums(is.na(interum_mat))==0]
  }
  return(list(class_rates = class_rates, best_feats = best_feats, best_feat_data = best_feat_data))
}

# Apply the function 
single_wtis <- check_tis(r2h, mr, ti60, ti61, ti62, ti63, ti64, rb, a, freqs, cells, n_trials, n_feats = 1)
single_wtis$best_feats # compared to Tree: 0.002 worse, LDA: 0.006 better, QDA: 0.001 better (didn't do RDA, too long)
pair_wtis <- check_tis(r2h, mr, ti60, ti61, ti62, ti63, ti64, rb, a, freqs, cells, n_trials, n_feats = 2)
pair_wtis$best_feats # compared to Tree: same, LDA: same, QDA: 0.001 better (didn't do RDA, too long)
trio_wtis <- check_tis(r2h, mr, ti60, ti61, ti62, ti63, ti64, rb, a, freqs, cells, n_trials, n_feats = 3)
trio_wtis$best_feats # compared to Tree: same, LDA: 0.001 better, QDA: same (didn't do RDA, too long)


################################## Check with just new features ###############################
check_tis_just_tis <- function(ti60, ti61, ti62, ti63, ti64, freqs, cells, n_trials, n_feats){
  # Construct single/pairs/trios of features
  feature.data <- data.frame(ti60, ti61, ti62, ti63, ti64) 
  meas_freq <- c(paste(rep(c("ti60", "ti61", "ti62", "ti63", "ti64"), each = length(freqs)), rep(freqs, 5)))
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
  class_rates <- data.frame(Features = all_feats, Tree = rep(0, mf_length), SE_Tree = rep(0, mf_length), LDA = rep(0, mf_length), SE_LDA = rep(0, mf_length), QDA = rep(0, mf_length), SE_QDA = rep(0, mf_length))
  cat("There are", mf_length, "features to go through.\n")
  
  for (k in 1:mf_length){
    if (k %% 100 == 0){
      cat("k:", k, "\n")
    }
    # Construct data frame to be used in classification algorithms
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
    tr <- rep(0, n_trials) # classification tree algorithm
    l <- rep(0, n_trials) # LDA
    q <- rep(0, n_trials) # QDA
    
    for (w in 1:n_trials){ # n_trials random splits of data
      set.seed(w) 
      
      # Get indicies of training set
      training_inds <- vector(length = length(cells))
      for (type in unique(cells)){
        inds <- which(cells == type)
        training_inds[sample(inds, 0.5*(length(inds)))] <- TRUE # randomly sample half of wells from each cell type
      }
      
      # Apply the three classification algorithms and out-of-sample predicitive accuracy
      tree.fit <- tree(cells ~ ., data = f.data, subset = training_inds)
      lda.fit <- lda(cells ~ ., data = f.data, subset = training_inds)
      qda.fit <- qda(cells ~ ., data = f.data, subset = training_inds)
      
      tree.pred <- predict(tree.fit, f.data[!training_inds, ],type = "class")
      tree.class <- tree.pred
      tr[w] <- mean(tree.class == f.data[!training_inds, 1])
      
      lda.pred <- predict(lda.fit, f.data[!training_inds, ])
      lda.class <- lda.pred$class
      l[w] <- mean(lda.class == f.data[!training_inds, 1])
      
      qda.pred <- predict(qda.fit, f.data[!training_inds, ])
      qda.class <- qda.pred$class
      q[w] <- mean(qda.class == f.data[!training_inds, 1])
      
    }
    class_rates[k, 2] <- mean(tr)
    class_rates[k, 3] <- sd(tr)/sqrt(n_trials)
    class_rates[k, 4] <- mean(l)
    class_rates[k, 5] <- sd(l)/sqrt(n_trials)
    class_rates[k, 6] <- mean(q)
    class_rates[k, 7] <- sd(q)/sqrt(n_trials)
  }
  
  # Tree in Column 2
  maxt <- class_rates[which.max(class_rates[, 2]), ]
  # LDA in Column 4
  maxl <- class_rates[which.max(class_rates[, 4]), ]
  # QDA in Column 6
  maxq <- class_rates[which.max(class_rates[, 6]), ]
  # Combine
  best_feats <- rbind(maxt, maxl, maxq)
  
  # Return data associated with the best features
  best_feat_data <- array(dim = c(nrow(ti60), n_feats, 3)) # for each classification method, matrix of size n_wells x n_feats
  for (i in 1:3){
    names_best_feats <- unlist(strsplit(as.character(best_feats[i, 1]), split = " "))
    rb_ind <- which(names_best_feats == "Rb")
    alpha_ind <- which(names_best_feats == "Alpha")
    interum_mat <- matrix(NA, nrow = nrow(ti60), ncol = length(names_best_feats)) # must accommodate Rb and Alpha which don't have frequencies, still want features in same order
    if (length(rb_ind) != 0){
      interum_mat[, rb_ind] <- feature.data[[which(meas_freq == "Rb")]]
    } 
    if (length(alpha_ind) != 0){
      interum_mat[, alpha_ind] <- feature.data[[which(meas_freq == "Alpha")]]
    }
    other_inds <- which(names_best_feats != "Rb"& names_best_feats != "Alpha")
    for (j in other_inds[seq(1, length(other_inds), by = 2)]){
      interum_mat[, j] <- feature.data[[which(meas_freq == paste(names_best_feats[j], names_best_feats[j+1]))]]
    }
    best_feat_data[,,i] <- interum_mat[, colSums(is.na(interum_mat))==0]
  }
  return(list(class_rates = class_rates, best_feats = best_feats, best_feat_data = best_feat_data))
}

# Apply the function 
single_wtis_just_tis <- check_tis_just_tis(ti60, ti61, ti62, ti63, ti64, freqs, cells, n_trials, n_feats = 1)
single_wtis_just_tis$best_feats # compared to Tree: .002 worse, LDA: .006 better, QDA: .001 better (didn't do RDA, too long)
pair_wtis_just_tis <- check_tis_just_tis(ti60, ti61, ti62, ti63, ti64, freqs, cells, n_trials, n_feats = 2)
pair_wtis_just_tis$best_feats # compared to Tree: .07 worse, LDA: .042 worse, QDA:  .001 better (didn't do RDA, too long)
trio_wtis_just_tis <- check_tis_just_tis(ti60, ti61, ti62, ti63, ti64, freqs, cells, n_trials, n_feats = 3)
trio_wtis_just_tis$best_feats # compared to Tree: .11 worse, LDA: .063 worse, QDA: .028 worse (didn't do RDA, too long)


