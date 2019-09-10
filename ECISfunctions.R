freqFeat <- function(freq.data){
  # Creates three different features based on time series data provided
  # Assumes one time index corresponds to twenty minutes
  #
  # Arguments:
  #    freq.data: data frame consisting of minimum of cell type, unique cell ids, times, frequencies, resistance values, Rb, and alpha values,
  #               order of data is important - assuming organized by cell type, well, frequency, and time (do have code to make sure this is case)
  # Returns: 
  #     List of three matrices (n_wells x n_freqs) corresponding to three different features, in order: 
  #         1. res2h - resistance at 2 hours 
  #         2. maxres - the maximum resistance across the time series 
  #         3. endofrun - average resistance over last 5 time points in time series 
  #     A unique value of 1. - 3. is obtained for each observation/well for each frequency of data collected
  #     This leads to matrices of dimension (n_wells x n_freqs)
  
  # Extract relevant dimensions
  n_freqs <- length(unique(freq.data$freq)) # number of unique frequencies
  n_ts <- length(unique(freq.data$time_index)) # length of time series (associated with each observation and frequency)
  n_ft <- n_freqs * n_ts # product of number of frequencies and length of time series
  n_wells <- nrow(freq.data)/n_ft # number of observations/wells of cells 
  
  # Create empty feature matrices
  # Columns correspond to all wells at particular frequency 
  res2h <- matrix(0, nrow = n_wells, ncol = n_freqs) 
  maxres <- matrix(0, nrow = n_wells, ncol = n_freqs)
  endofrun <- matrix(0, nrow = n_wells, ncol = n_freqs)
  
  # Make sure data is ordered correctly, assuming organized by cell type, well, frequency, and time
  freq.data.ordered <- freq.data[order(freq.data$cell, freq.data$well, freq.data$freq, freq.data$time_index), ]
  
  # Generate the features
  for (w in 1:n_wells){ # For each well,
    start <- n_ft*(w-1) + 1 # start in correct row of data frame (associated with first instance of well w)
    fin <- start + (n_ft - 1) # end in correct row of data frame
    X <- matrix(freq.data.ordered[start:fin, "res"], n_ts, n_freqs) # only need resistance measurements for feature generation, time series of length n_ts, each column for each frequency 
    for (f in 1:n_freqs){ # For each frequency,
      res2h[w, f] <- rollmean(X[, f], k = 5, fill = FALSE)[6] # row/time index 6 corresponds to 2 hour mark, column indicates frequency
      maxres[w, f] <- max(rollmean(X[, f], k = 5, fill = FALSE)) # use max of SMA time series instead of original
      endofrun[w, f] <- mean(X[(n_ts - 4):n_ts, f]) 
    }
  }
  return(list(res2h, maxres, endofrun))
}

nbestFeat <- function(r2h, mr, eor, rb, a, freqs, cells, n_trials, n_feats){
  # Uses three different classification techniques on data/features provided to find
  # best "n" features for cell line classification (for each technique).
  # Currently, only allows n = 1, 2, 3.
  # Calculates out-of-sample predictive accuracy as measure of goodness (of feature).
  # N best features selected as majority vote over n_trials different random splits of the data.
  # 
  # Arguments:
  #     r2h: matrix of res2h values (n_wells x n_freqs)
  #     mr: matrix of maxres values (n_wells x n_freqs)
  #     eor: matrix of endofrun values (n_wells x n_freqs)
  #     rb: vector of Rb values, same across all frequencies (n_wells x 1)
  #     a: vector of alphda values, same across all frequencies (n_wells x 1)
  #     freqs: vector of unique frequencies (n_freqs x 1)
  #     cells: vector of cell types associated with each well, same across all frequencies (n_wells x 1)
  #     n_trials: number of trials (random splits of data to perform) 
  #     n_feats: number of best features to consider - 
  #               n = 1: single best
  #               n = 2: pair of best
  #               n = 3: trio of best
  # Returns:
  #     List of 3 objects:
  #         1. class_rates: full matrix of features and associated out-of-sample predicitive accuracy rates and s.e. for each classification technique
  #         2. best_feats: (3 x n_feats) matrix of best classification feature according to row 1: classification tree, row 2: LDA, row 3: QDA
  #         3. best_feat_data: (n_wells x n_feats x 3) array of the data associated with the best features for each classification method
  #                                 [,,1]: classification trees
  #                                 [,,2]: LDA
  #                                 [,,3]: QDA
  
  # Construct single/pairs/trios of features
  feature.data <- data.frame(r2h, mr, eor, rb, a) # r2h, mr, and eor matrices with columns corresponding to different frequencies
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
  class_rates <- data.frame(Features = all_feats, Tree = rep(0, mf_length), SE_Tree = rep(0, mf_length), LDA = rep(0, mf_length), SE_LDA = rep(0, mf_length), QDA = rep(0, mf_length), SE_QDA = rep(0, mf_length))
  
  for (k in 1:mf_length){
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
      qda.fit <- try(qda(cells ~ ., data = f.data, subset = training_inds), silent = TRUE) # QDA

      tree.pred <- predict(tree.fit, f.data[!training_inds, ],type = "class")
      tree.class <- tree.pred
      tr[w] <- mean(tree.class == f.data[!training_inds, 1])

      lda.pred <- predict(lda.fit, f.data[!training_inds, ])
      lda.class <- lda.pred$class
      l[w] <- mean(lda.class == f.data[!training_inds, 1])

      if (inherits(qda.fit, "try-error")){ # Sometimes collinearity - don't use that feature combo
        #cat(all_feats[k], "\n")
        q[w] <- 0
      } else{
        qda.pred <- predict(qda.fit, f.data[!training_inds, ])
        qda.class <- qda.pred$class
        q[w] <- mean(qda.class == f.data[!training_inds, 1])
      }
      
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
  possible_feats <- c("res2h", "maxres","endofrun", "Rb", "Alpha")
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

rdaScan <- function(r2h, mr, eor, rb, a, freqs, cells, n_trials, n_feats, rhos){
  # Performs scan of fine scale of rho_1 values to find best one for RDA analysis
  # Returns array of RDA out-of-sample predictive accuracy rates and standard errors
  # for each n_feats combinations of features
  #
  # Arguments:
  #     r2h: matrix of res2h values (n_wells x n_freqs)
  #     mr: matrix of maxres values (n_wells x n_freqs)
  #     eor: matrix of endofrun values (n_wells x n_freqs)
  #     rb: vector of Rb values, same across all frequencies (n_wells x 1)
  #     a: vector of alphda values, same across all frequencies (n_wells x 1)
  #     freqs: vector of unique frequencies (n_freqs x 1)
  #     cells: vector of cell types associated with each well, same across all frequencies (n_wells x 1)
  #     n_trials: number of trials (random splits of data to perform) 
  #     n_feats: number of features (pairs = 2, trios = 3) search over to find best rho_1
  #     rhos: sequence of rho values to test
  # Returns:
  #    class_rates_RDA: array of RDA out-of-sample predictive accuracy rates and standard errors for each 
  #                     n_feats combination of features and each rho_1 values
  
  # Collect feature data
  feature.data <- data.frame(r2h, mr, eor, rb, a) # r2h, mr, and eor matrices with columns corresponding to different frequencies
  meas_freq <- c(paste(rep(c("res2h", "maxres","endofrun"), each = length(freqs)), rep(freqs, 3)), "Rb", "Alpha")
  # Specifcy range of rho_1 values want to test
  rho_1 <- rhos
  # Construct feature combinations for testing rho_1 values
  if (n_feats == 2){
    feat_inds <- t(combn(1:length(meas_freq),2))
    mf_length <- nrow(feat_inds)
    all_feats <- paste(meas_freq[feat_inds[, 1]], meas_freq[feat_inds[, 2]])
  } else {
    feat_inds <- t(combn(1:length(meas_freq),3))
    mf_length <- nrow(feat_inds)
    all_feats <- paste(meas_freq[feat_inds[, 1]], meas_freq[feat_inds[, 2]], meas_freq[feat_inds[, 3]])
  }
  class_rates_RDA <- array(dim = c(mf_length, 3, length(rho_1)))
  colnames(class_rates_RDA) <- c("Features", "RDA", "SE_RDA")
  class_rates_RDA[,1,] <- all_feats
 
  # Do the testing
  list_results <- mclapply(rho_1, function(x){
    small_scale <- matrix(nrow = mf_length, ncol = 2)
    for (k in 1:mf_length){
      feature1 <- feature.data[, feat_inds[k, 1]]
      feature2 <- feature.data[, feat_inds[k, 2]]
      if (n_feats == 3){
        feature3 <- feature.data[, feat_inds[k, 3]]
        f.data <- data.frame(cells, feature1, feature2, feature3)
      } else{
        f.data <- data.frame(cells, feature1, feature2)
      }
      
      # Create empty matrix to store out-of-sample predictive accuracy for each rho_1
      rda <- vector(length = n_trials) 
      
      for (w in 1:n_trials){ # n_trials random splits of data
        set.seed(w)
        
        # Get indicies of training set
        training_inds <- vector(length = length(cells))
        for (type in unique(cells)){
          inds <- which(cells == type)
          training_inds[sample(inds, 0.5*(length(inds)))] <- TRUE # randomly sample half of wells from each cell type
        }  
        
        trainY <- f.data[training_inds, -1]
        testY <- f.data[!training_inds, -1]
        
        res.RDA <- gels_RDA(trainY, as.numeric(factor(cells[training_inds])), testY, estim="class", r1 = x)
        rda[w] <- sum(res.RDA$predict==as.numeric(factor(cells[!training_inds])))/sum(!training_inds) # out-of-sample predictive accuracy
      }
      small_scale[k, 1] <- mean(rda)
      small_scale[k, 2] <- sd(rda)/sqrt(n_trials)
    }
    return(small_scale)
  }, mc.cores = detectCores())

  for (i in 1:length(rho_1)){
    class_rates_RDA[, 2:3, i] <- list_results[[i]] # mean and s.e. of out-of-sample predictive accuracies for RDA with particular rho_1
  }

  # What are the best features associated with the best overall rho_1? (where do we find the highest out-of-sample predictive accuracy?)
  best_inds <- which(class_rates_RDA[, 2, ] == max(class_rates_RDA[, 2, ]), arr.ind = TRUE) 
  if (length(dim(best_inds)) == 0){ # if only one best
    best_feats <- c(class_rates_RDA[best_inds[1], , best_inds[2]], sprintf("%.2f", rho_1[best_inds[2]]))
  } else{ # if more than one best, just chose first
    best_feats <- c(class_rates_RDA[best_inds[1, 1], , best_inds[1, 2]], sprintf("%.2f", rho_1[best_inds[1, 2]]))
  }
  
  # Return data associated with the best feature
  #possible_feats <- c("res2h", "maxres","endofrun", "Rb", "Alpha")
  names_best_feats <- unlist(strsplit(as.character(best_feats[1]), split = " "))
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
  best_feat_data <- interum_mat[, colSums(is.na(interum_mat))==0]

  return(list(class_rates_RDA = class_rates_RDA, best_feats = best_feats, best_feat_data = best_feat_data))
}

gels_RDA<-function(trainY, groupings, testY,
                   estim="class", r1, r2 = 0,
                   Maxmin.ratio=0.01, Maxmin.ratio2=0.01, nsteps=5,nsteps2=5){
  
  # Function from Ines Wilms and Stephanie Aerts - slightly modified so that rho_2 = 0 fixed 
  # Removed code associated with rho_2
  
  ######################
  # INPUTS
  ######################
  # trainY : Data matrix of dimension N X p (N: total sample size, p: number of variables)
  # groupings : numeric vector of dimension N containing the group membership of the observations
  # testY : Data matrix of dimension m X p  containing the m new observations to classify
  # estim : character string for the choice of location and covariance estimates
  #         either "class" for the sample estimators, or "pairwCorrK" for the cordinatewise median and the pairwise covariance matrices
  #         based on Kendall's correlation
  # Maxmin.ratio : ratio between the lower bound and the upper bound of the grid of values for the parameter rho1
  # Maxmin.ratio2 : ratio between the lower bound and the upper bound of the grid of values for the parameter rho2
  # nsteps: number of parameter values in the grid search for the optimal rho1
  # nsteps2: number of parameter values in the grid search for the optimal rho2
  
  ######################
  # OUTPUTS
  ######################
  # estim: choice of location and covariance estimates 
  # mean_list : list of mean estimates for each group
  # S.list : list of covariance matrix estimates for each group
  # S.RDA: list of regularized covariance matrix estimates for each group
  # rho1_max: maximum value of the grid of values for rho1
  # rho1_opt : optimal regularization parameter lambda1 minimizing the BIC
  # rho1_max: maximum value of the grid of values for rho2
  # rho1_opt : optimal regularization parameter lambda2 minimizing the BIC
  # classes: predicted groups for the observation of testY
  
  ######################
  # LIBRARIES
  ######################
  # library(huge)
  # library(pcaPP)
  # library(robustbase)
  
  ######################
  # AUXILIARY FUNCTIONS
  ######################
  # CellwiseCov.R
  # BICfct.R
  # predictDA.R
  
  ######################
  # ADDITIONAL FUNCTIONS
  ######################
  is.invertible <- function(m) class(try(solve(m),silent=T))=="matrix"
  
  RDAest<-function(S.list,S.pooled, rho_1, rho_2,n,p){
    temp1<-lapply(S.list,function(x, rho_1, S.pooled)(1-as.numeric(rho_1))*x + as.numeric(rho_1)*as.matrix(S.pooled), rho_1=rho_1, S.pooled=S.pooled)
    sigmaRDA<- lapply(temp1, function(x, rho_2,p)(1-rho_2)*x +rho_2/p *sum(diag(x))* diag(p),rho_2=rho_2, p=p)
    return(list(sigmaRDA=sigmaRDA, rho_1=rho_1, rho_2=rho_2))
  }
  
  BIC.RDA <-function( rho1,rho2, S.list,S.pooled, n, p){
    FIT<-RDAest(S.list=S.list, S.pooled, n=n,p=p,rho_1=rho1,rho_2=rho2) 
    if(Reduce('+',lapply(FIT$sigmaRDA,  is.invertible ))==length(S.list)){
      BIC_value<- BICfct(A=S.list, B= lapply(FIT$sigmaRDA,solve), n=n,p=p)
    }else{
      BIC_value=NA
    }
    return(BIC_value)
  }
  
  
  ##########
  # START
  ##########
  trainY<-split(data.frame(trainY),groupings)
  K<-length(unique(groupings))
  n<-unlist(lapply(trainY,nrow))
  p<-dim(trainY[[1]])[2]
  
  # Estimation of the means and covariance matrices
  if(estim=="class"){
    mean_list <- lapply(trainY, function(u)apply(u,2,mean))
    S.list<-lapply(trainY, cov)
  }else if(estim=="pairwCorrK"){
    mean_list<- lapply(trainY,function(t) apply(t,2,median))
    S.list<- lapply(trainY, kendall.transformed)
  }
  S.pooled<- Reduce('+',mapply(function(s,n)n*s,s=S.list,n=n,SIMPLIFY=F))/(sum(n)-K)
  
  # Grid of regularization parameters
  if (r1 == 99){
    rho1_max<- 1 
    rho1_seq<-sort(exp(seq(log(rho1_max),log(Maxmin.ratio*rho1_max),length.out=nsteps)),decreasing=T)
    r1<-list(rho1_l=sort(rep(rho1_seq),decreasing=T), rho2_l=rep(rho2,nsteps))
    
    # Selection of regularization parameters
    morearg.plug=list(S.list=S.list, S.pooled=S.pooled,n=n,p=p)
    BICval.RDA<-mapply(BIC.RDA,r1$rho1_l, r1$rho2_l, MoreArgs = morearg.plug,SIMPLIFY=F)
    rho1_opt.BIC<- r1$rho1_l[which.min(BICval.RDA)]
    rho1 <- rho1_opt.BIC
  } else{
    rho1 <- r1
    rho2 <- r2
  }
  # Final solution with optimal parameters
  RDA.BIC<- RDAest(S.list=S.list,S.pooled= S.pooled,
                   rho_1=rho1, rho_2=r2,
                   n=n,p=p)
  classes<-predict.DA(testY, means=mean_list,theta=lapply(RDA.BIC$sigmaRDA, solve),type="QDA",priors=n/sum(n))$cl
  
  return(list( estim=estim,mean_list=mean_list,S.list=S.list, 
               S.RDA= RDA.BIC$sigmaRDA, 
               rho1_opt=rho1,
               rho2_opt=rho2,
               predict = classes))
}


predict.DA<- function (means, theta, priors, x, type) 
{
  # Function from Ines Wilms and Stephanie Aerts - No modifications mode.
  
  ######################## 
  # INPUTS
  ########################
  # means: list of K mean estimates
  # theta : either a list of K precision matrices (when QDA) or one single common precision matrix estimate
  # priors : vector of prior probabilities of each group
  # x : new data for which the groups are to predict
  # type: either QDA or LDA
  ######################## 
  # OUTPUTs
  ########################
  # class: class prediction for each observation
  # posterior: discriminant values
  
  if(type=="QDA"){
    probs <- apply(x, 1,function(x) mapply( function(mean, omega,priors) {0.5* log(det(as.matrix(omega)))  -
        0.5* as.double(unlist(x-mean)%*% as.matrix(omega) %*% unlist(x-mean)) +
        log(priors)}, omega=theta, mean=means,priors=priors))
  }
  if(type=="LDA"){
    probs <- apply(x, 1,function(x) mapply( function(mean,priors) {0.5* log(det(as.matrix(theta)))  -
        0.5* as.double(unlist(x-mean)%*% as.matrix(theta) %*% unlist(x-mean)) +
        log(priors)},mean=means,priors=priors))
  }
  cl = apply(probs, 2,function(p)which.max(p))
  score= unlist(lapply(1:nrow(x), function(i) probs[cl[i],i]))
  return(list(cl=cl,score=score))
}

getExample <- function(cells, data, class_type, rho_1 = NA, train){
  # Function for generating an example of classification method fit
  # Needed for various figures in manuscript
  #
  # Arguments:
  #     cells: vector of cell types associated with each well, same across all frequencies (n_wells x 1)
  #     data:  (n_wells x 2) matrix of data to use in plotting
  #     class_type:  classification method, supports: ("tree", "lda", "qda", "rda")
  #     rho_1: rho_1 value to use with "rda", defualt is NA 
  #     train:  logical, train the contours on subset of data (TRUE) or use entire data (FALSE)
  # Returns:
  #     List of four vector values needed to construct plots (nd.x, nd.y, prd, pch_vals)
  
  set.seed(2) 
  
  # Construct data
  colnames(data) <- c("X1", "X2")
  f.data <- data.frame(cells, data)
  n <- length(cells)
  
  if (train){ # If using training set...
    # Get indicies of training set
    training_inds <- vector(length = n)
    for (type in unique(cells)){
      inds <- which(cells == type)
      training_inds[sample(inds, 0.5*(length(inds)))] <- TRUE # randomly sample half of wells from each cell type
      pch_vals <- vector(length = n)
      pch_vals[training_inds] <- 1 # training data have 'o' symbol
      pch_vals[!training_inds] <- 4 # other data points have 'x' symbol, not used to construct contours
    }
  } else{ # Otherwise...
    training_inds <- 1:n # use entire data set
    pch_vals <- rep(1, n) # all data have same 'o' symbol
  }
  
  # Set up grid of values
  np <- 300
  nd.x <- seq(from = min(data[, 1]), to = max(data[, 1]), length.out = np)
  nd.y <- seq(from = min(data[, 2]), to = max(data[, 2]), length.out = np)
  nd <- expand.grid(X1 = nd.x, X2 = nd.y)
  
  # Perform fit and prediction
  if (class_type == "tree"){
    fit <- tree(cells ~ ., data = f.data, subset = training_inds)
  } else if (class_type == "lda"){
    fit <- lda(cells ~ ., data = f.data, subset = training_inds)
  } else if (class_type == "qda"){
    fit <- qda(cells ~ ., data = f.data, subset = training_inds)
  } else if (class_type == "rda"){
    trainY <- f.data[training_inds, -1]
    fit <- gels_RDA(trainY, as.numeric(factor(cells[training_inds])), nd, estim="class", r1 = rho_1)
  }
  
  if (class_type == "tree"){
    prd <- as.numeric(predict(fit, newdata = nd, type = "class"))
  } else if (class_type == "lda" | class_type == "qda"){
    prd <- as.numeric(predict(fit, newdata = nd)$class)
  } else if (class_type == "rda"){
    prd <- as.numeric(fit$predict)
  }
  
  return(list(nd.x = nd.x, nd.y = nd.y, prd = prd, pch_vals = pch_vals))
}

createLab <- function(best_feats){
  # Convert strings indicating best features to nice format for figures and tables
  #
  # Arguments:
  #    best_feats: vector of strings indicating best features in data set 
  # Returns:
  #    feat_labs: vector of strings indicating best features in data set in nice format
  
  names_best_feats <- unlist(strsplit(as.character(best_feats), split = " ")) 
  rb_alpha_inds <- which(names_best_feats == "Rb"| names_best_feats == "Alpha")
  feat_labs <- rep(NA, length(names_best_feats)) # must accommodate Rb and Alpha which don't have frequencies
  if (length(rb_alpha_inds) != 0){
    feat_labs[rb_alpha_inds] <- names_best_feats[rb_alpha_inds]
  }
  other_inds <- which(names_best_feats != "Rb"& names_best_feats != "Alpha")
  for (i in other_inds[seq(1, length(other_inds), by = 2)]){
    feat_labs[i] <- createLab_ind(names_best_feats[i], names_best_feats[i+1])
  }
  feat_labs <- feat_labs[!is.na(feat_labs)]
  return(feat_labs)
}

createLab_ind <- function(feat_name, freq_name){
  # Called internally from "createLab"
  # Convert strings indicating best features to nice format for figures and tables (for one)
  # Assuming not passing "Rb" or "Alpha"
  #
  # Arguments:
  #    feat_name: string indicating name of best feature
  #    freq_name: string indicating frequency associated with feature
  # Returns:
  #   feat_lab: string used for axis label, represents best feature in nice format
  
  if (feat_name == "res2h"){
    lab_name <- "R2h"
  } else if (feat_name == "maxres"){
    lab_name <- "MR"
  } else if (feat_name == "endofrun"){
    lab_name <- "EOR"
  } else if (grepl("ti", feat_name)){
    lab_name <- paste("TI", substr(feat_name, 3, 4), sep = "")
  }
  feat_lab <- paste(lab_name, "@", freq_name, "Hz")

  return(feat_lab = feat_lab)
}

mainAnalysis <- function(freq.data){
  # Executes analysis performed in IJBS paper
  # Constructs tables and figures in IJBS paper
  #
  # Arguments:
  #    freq.data: data frame consisting of minimum of cell type, unique cell ids, times, frequencies, resistance values, Rb, and alpha values,
  #               order of data is important - assuming organized by cell type, well, frequency, and time (do have code to make sure this is case)
  # Returns:
  #    List of analysis output
  
  # Construct the frequency features
  freqfeat <- freqFeat(freq.data)
  res2h <- freqfeat[[1]]
  maxres <- freqfeat[[2]]
  endofrun <- freqfeat[[3]]
  
  # Extract key information from data - useful in analysis (see below)
  n_freqs <- length(unique(freq.data$freq)) # number of unique frequencies
  n_ts <- length(unique(freq.data$time_index)) # length of time series (associated with each observation and frequency)
  n_ft <- n_freqs * n_ts # product of number of frequencies and length of time series
  inds <- seq(1, nrow(freq.data), by = n_ft) # need one index associated with each well since Rb, alpha and cell type same across times and frequencies
  rb <- freq.data$rb[inds]
  alpha <- freq.data$alpha[inds]
  cells <- freq.data$cell[inds]
  freqs <- unique(freq.data$freq)
  n_trials <- 20 
  
  # Find single best feature (7-8 sec. for gel, )
  singlefeat <- nbestFeat(res2h, maxres, endofrun, rb, alpha, freqs, cells, n_trials, n_feats = 1)
  cat("The single best feature has been found.\n")
  
  # Find pair of best features (2 min. for gel, )
  pairfeat <- nbestFeat(res2h, maxres, endofrun, rb, alpha, freqs, cells, n_trials, n_feats = 2)
  cat("The best pair of features has been found.\n")
  
  # Find trio of best features (17 min. for gel, )
  triofeat <- nbestFeat(res2h, maxres, endofrun, rb, alpha, freqs, cells, n_trials, n_feats = 3)
  cat("The best trio of features has been found.\n")
  
  # Perform RDA scan, assess best rho_1 value (and features associated with that)...
  # For pairs of features...(2 h. for gel, )
  rhos <- seq(0.05, 0.95, by = 0.05)
  pair_rdascan <- rdaScan(res2h, maxres, endofrun, rb, alpha, freqs, cells, n_trials, n_feats = 2, rhos)
  cat("The RDA scan on all pairs of features has been completed.\n")
  
  # And trios of features
  # Same as above: rhos <- seq(0.05, 0.95, by = 0.05)
  trio_rdascan <- rdaScan(res2h, maxres, endofrun, rb, alpha, freqs, cells, n_trials, n_feats = 3, rhos)
  cat("The RDA scan on all trios of features has been completed.\n")
  
  return(list(cells = cells, rhos = rhos, singlefeat = singlefeat, pairfeat = pairfeat, triofeat = triofeat, pair_rdascan = pair_rdascan, trio_rdascan = trio_rdascan))
  
}


