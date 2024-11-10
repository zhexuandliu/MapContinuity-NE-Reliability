########################################
# Nov 1, 2024, Zhexuan Liu #############
# Verify LOO assumption ################
########################################

library(Rfast)
library(Rtsne)
library(uwot)

generate_data = function(n, d = 2, data_type){
  if (data_type == '2GMM'){
    dist =  6 # set component distance
    dim = d # set data dimension
    X = matrix(rnorm(dim * n), nrow = n)
    X[1:round(n/2),1] = X[1:round(n/2),1] + dist/2
    X[(round(n/2)+1):n,1] = X[(round(n/2)+1):n,1] - dist/2
    return(X)
  }else if (data_type == 'swissroll'){
    dataset_name = 'swissroll'
    generate_swiss_roll <- function(n_samples = 1500, noise = 0.05) {
      t = runif(n_samples, 0, 3*pi)
      x = t * cos(t)
      y = t * sin(t)
      z = runif(n_samples, -1, 1) * 4
      # Add noise
      x = x + rnorm(n_samples, mean = 0, sd = noise)
      y = y + rnorm(n_samples, mean = 0, sd = noise)
      z = z + rnorm(n_samples, mean = 0, sd = noise)
      list(x = x, y = y, z = z, t = t)
    }
    X = generate_swiss_roll(n_sample = n)
    return(cbind(X$x, X$y, X$z))
  }else if (data_type == 'panc8'){
    load('./data/VerifyLOO/datasets/panc8_1&2.RData')
    X = t(data.X1)
    X = X[sample(c(1:dim(X)[1]), n, replace = FALSE),]
  }else if (data_type == 'ifnb'){
    load('./datasets/ifnb.RData')
    X = t(data.X1)
    X = X[sample(c(1:dim(X)[1]), n, replace = FALSE),]
  }else if (data_type == 'brain'){
    load('./datasets/brain.RData')
    X = t(data.X1)
    X = X[sample(c(1:dim(X)[1]), n, replace = FALSE),]
  }else if (data_type == 'cifar10'){
    X = read.csv("./datasets/cifar10.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE, na.strings = "NA") # data can be downloaded from https://drive.google.com/drive/folders/1mrV0ROvJ803k_3B7FVqNM5u6HaLd8FPv?usp=share_link
    X = as.matrix(X)
    X = X[sample(c(1:dim(X)[1]), n, replace = FALSE),]
  }else{
    print('data type not supported')
  }
}

loo_diff = function(X, perplexity, method = 'tsne'){
  # generate embedding based on original data with one random point added
  if ((dim(X)[1] - 2) < 3 * perplexity){
    return(NA)
  }
  if (method == 'umap'){
    Y = umap(X, n_neighbors = perplexity, n_epoch = 2000)
  }else{
    Y = Rtsne(X, perplexity = perplexity, theta = 0.0, max_iter = 3000)
    Y = Y$Y
  }

  # embedding before adding one point
  X = X[-1, ]
  Y_tilde = Y[-1,]
  if (method == 'umap'){
    initial_embd = umap(X, n_neighbors = perplexity, init = Y_tilde)
  }else{
    initial_embd = Rtsne(X, perplexity = perplexity, theta = 0.0, Y_init = Y_tilde)
    initial_embd = initial_embd$Y
  }
  diff = Y_tilde - initial_embd
  epsilon = 1 / (sum(initial_embd ^ 2))**0.5 * (sum(diff ^ 2))**0.5
  return(epsilon)
}

# We verified the LOO assumption for various datasets and different perplexities. The results are stored under ./data/VerifyLOO/epsilon.
