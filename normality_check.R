load("gel_data.Rdata")
load("bsa_data.Rdata")
source("ECISfunctions.R")

# Construct the frequency features
freq.data <- gel.data # bsa.data
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

# Construct histograms to assess normality within each cell type for each feature
par(mfrow = c(3, 5))
for (j in 1:ncol(maxres)){
  for (i in unique(cells)){
    hist(maxres[cells == i, j], xlab = paste(i), main = paste("MR Column:", j))
  }
}
# Okay for maxres

for (j in 1:ncol(res2h)){
  for (i in unique(cells)){
    hist(res2h[cells == i, j], xlab = paste(i), main = paste("R2h Column:", j))
  }
}
# Okay for res2h

for (j in 1:ncol(endofrun)){
  for (i in unique(cells)){
    hist(endofrun[cells == i, j], xlab = paste(i), main = paste("EOR Column:", j))
  }
}
# Okay for endofrun

for (i in unique(cells)){
  hist(alpha[cells == i], xlab = paste(i), main = paste("Alpha"))
}
# Okay for alpha

for (i in unique(cells)){
  hist(rb[cells == i], xlab = paste(i), main = paste("Rb"))
}
# Okay for rb

## OVERALL: Normality within class reasonable assumption for all features - some cases where dat might be slightly left-skewed or right-skewed but nothing extreme. 