###############################################################################
############################## Singularity Score ##############################
###############################################################################

#' gradient_compute
#'
#' computes gradients for t-SNE loss function
#'
#' @param Y Matrix. t-SNE embedding.
#' @param P Matrix. Similarity matrix in t-SNE.
#' @return Gradients for t-SNE loss function at embedding Y.
#' @examples
#' gradient_compute(Y, P)
#' @export
#' @import Rfast
gradient_compute = function(Y, P){
  Y = c(t(Y))
  n = dim(P)[1]
  n2 = dim(P)[2]
  n_Y = length(Y)
  if ((n!=n2) | (2*n!=n_Y)){
    stop('Input dimensions do not match!')
  }

  YMatmatrix = matrix(rep(c(matrix(matrix(Y, n, 2, byrow = TRUE), 1, 2*n)), n), 2*n, n, byrow = TRUE)
  YMat = matrix(Y, 2, n)
  YnormSq = matrix(colSums(YMat ** 2), n, 1)
  YDiffDistSq = YnormSq %*% matrix(1,1,n) + matrix(1,n,1) %*% t(YnormSq) - 2 * t(YMat) %*% YMat
  YDiffDistSq_double = matrix(rep(c(matrix(YDiffDistSq)), each = 2), 2*n, n)
  P_double = matrix(rep(c(matrix(P)), each = 2), 2*n, n)
  deno = sum(1 / (1 + YDiffDistSq)) - n

  I1 = 4 * P_double * (matrix(rep(Y,n),2*n,n) - YMatmatrix) / (1 + YDiffDistSq_double)
  I2 = (-1) * 4 * (matrix(rep(Y,n),2*n,n) - YMatmatrix) / deno / ((1 + YDiffDistSq_double) ^ 2) * sum(P)

  G = rowSums(I1 + I2)

  return(G)
}

#' hessian_compute
#'
#' computes the Hessian matrix for t-SNE loss function
#'
#' @param Y Matrix. t-SNE embedding.
#' @param P Matrix. Similarity matrix in t-SNE.
#' @return Hessian matrix for t-SNE loss function at embedding Y.
#' @examples
#' hessian_compute(Y, P)
#' @export
#' @import Rfast
hessian_compute = function(Y, P){
  Y = c(t(Y))
  n = dim(P)[1]
  n2 = dim(P)[2]
  n_Y = length(Y)
  if ((n!=n2) | (2*n!=n_Y)){
    stop('Input dimensions do not match!')
  }

  blockGramMatmatrix = matrix(nrow=2*n, ncol=2*n)
  blockGramMatmatrix2 = matrix(1, nrow=2*n, ncol=2*n)
  indices1 = seq(from=1, to=2*n, by=2)
  indices2 = seq(from=2, to=2*n, by=2)
  blockGramMatmatrix2[indices1, indices2] = 0
  blockGramMatmatrix2[indices2, indices1] = 0
  I2 = matrix(nrow=2*n, ncol=2*n)
  I4 = matrix(nrow=2*n, ncol=2*n)
  I5 = matrix(nrow=2*n, ncol=2*n)

  YMat = matrix(Y, 2, n)
  YnormSq = matrix(colSums(YMat ** 2), n, 1)
  YDiffDistSq = YnormSq %*% matrix(1,1,n) + matrix(1,n,1) %*% t(YnormSq) - 2 * t(YMat) %*% YMat
  deno = sum(1 / (1 + YDiffDistSq)) - n
  YMatmatrix = matrix(rep(c(matrix(matrix(Y, n, 2, byrow = TRUE), 1, 2*n)), n), 2*n, n, byrow = TRUE)
  YDiffDistSq_double = matrix(rep(c(matrix(YDiffDistSq)), each = 2), 2*n, n)

  I1 = ((-4) * P / (1 + YDiffDistSq)) %x% diag(2)
  I3 = - 16 * sum(P) / deno**2 * matrix(c(rowSums((matrix(rep(Y,n),2*n,n) - YMatmatrix) / ((1 + YDiffDistSq_double) ^ 2))), 2*n, 1) %*% matrix(c(rowSums((matrix(rep(Y,n),2*n,n) - YMatmatrix) / ((1 + YDiffDistSq_double) ^ 2))), 1, 2*n)
  I4 = (4 * sum(P) / deno / ((1 + YDiffDistSq)**2)) %x% diag(2)

  for (coordinate_idx1 in c(1:2)){
    for (coordinate_idx2 in c(1:2)){
      u = YMat[coordinate_idx1,]
      v = YMat[coordinate_idx2,]
      indices1 = seq(from=coordinate_idx1, to=2*n, by=2)
      indices2 = seq(from=coordinate_idx2, to=2*n, by=2)
      tmp = matrix(u,n,1) %*% matrix(v,1,n)
      blockGramMatmatrix[indices1, indices2] = matrix(rep(diag(tmp),n),n,n) + matrix(rep(diag(tmp),each=n),n,n) - tmp - t(tmp)
      I2[indices1, indices2] = (8 * P) / ((1 + YDiffDistSq)**2) * blockGramMatmatrix[indices1, indices2]
      I5[indices1, indices2] = (-1) * (16 * sum(P)) / deno / ((1 + YDiffDistSq)**3) * blockGramMatmatrix[indices1, indices2]
    }
  }

  H = I1 + I2 + I3 + I4 + I5
  for (coordinate_idx1 in c(1:2)){
    for (coordinate_idx2 in c(1:2)){
      indices1 = seq(from=coordinate_idx1, to=2*n, by=2)
      indices2 = seq(from=coordinate_idx2, to=2*n, by=2)
      diag(H[indices1,indices2])=0
    }
  }
  tmp = apply(array(H, dim=c(2,n,2,n)), c(1,2,3), sum)
  for (i in c(1:n)){
    H[c(2*i-1, 2*i), c(2*i-1, 2*i)] = (-1) * tmp[,i,]
  }

  return (H)
}

#' singularity_score_compute
#'
#' computes the singularity score to diagnose fracture-inducing discontinuity
#'
#' @param Y Matrix. t-SNE embedding.
#' @param P Matrix. Similarity matrix in t-SNE.
#' @return Hessian matrix for t-SNE loss function at embedding Y.
#' @examples
#' singularity_score_compute(Y, P)
#' @export
#' @import Rfast
singularity_score_compute = function(Y, P){
  eigen_score_pointwise = function(Hessian, ind){
    h = Hessian[c(2*ind-1, 2*ind), c(2*ind-1, 2*ind)]
    return(1/min(eigen(h)$values))
  }
  H = hessian_compute(Y, P)
  s_score = sapply(c(1:dim(Y)[1]), function(x) eigen_score_pointwise(Hessian = H, x))
  return(s_score)
}


###############################################################################
############################# Perturbation Score ##############################
###############################################################################

#' LOO_loss_fast
#'
#' get t-SNE LOO loss for specified point
#'
#' @param yy Vector. Embedding of the point that we want to calculate loss with.
#' @param Y Matrix. Embedding matrix of other points.
#' @param YDistSqP1 Matrix. Squared distance matrix plus 1 of 'rbind(c(0,0), Y)' for fast computation.
#' @param PMat Matrix. Simalarity matrix (P matrix) of all points (with the point we want to calculate loss with as the first point).
#' @return The LOO loss for the point yy.
#' @examples
#' LOO_loss_fast(yy, Y, YDistSqP1, Mat)
#' @export
#' @import Rfast
LOO_loss_fast = function(yy, Y, YDistSqP1, PMat){
  n = dim(Y)[1] + 1
  yy_Y_dist_sqP1 = 1 + rowsums((Y - matrix(yy, nrow = n-1, ncol = 2, byrow = TRUE))**2)
  YDistSqP1[1,2:n] = yy_Y_dist_sqP1
  YDistSqP1[2:n,1] = yy_Y_dist_sqP1
  I1 = 2 * (PMat[1,] %*% log(YDistSqP1[,1]))[1,1]
  I2 = log(sum(rowsums(1/YDistSqP1))-n)
  return(I1 + I2)
}

#' LOO_loss_compute
#'
#' get t-SNE LOO loss for specified point
#'
#' @param yy Vector. Proposed embedding point of the point x that we want to calculate loss with.
#' @param Y Matrix. Embedding matrix of X.
#' @param PMat Matrix. Simalarity matrix (P matrix) of 'rbind(x,X)'.
#' @return The LOO loss for the point yy.
#' @examples
#' LOO_loss_compute(yy, Y, Mat)
#' @export
#' @import Rfast
LOO_loss_compute = function(yy, Y, PMat){
  n = dim(Y)[1] + 1
  YDistSqP1 = as.matrix(Dist(rbind(yy, Y), square = TRUE)) + 1
  I1 = 2 * (PMat[1,] %*% log(YDistSqP1[,1]))[1,1]
  I2 = log(sum(rowsums(1/YDistSqP1))-n)
  return(I1 + I2)
}

#' gradient_compute_pointwise
#'
#' computes gradients for t-SNE loss function pointwise
#'
#' @param yy Vector. Embedding of the point that we want to calculate loss with.
#' @param Y Matrix. Embedding matrix of other points.
#' @param P Matrix. Similarity matrix in t-SNE.
#' @return Gradients for t-SNE loss function at embedding yy.
#' @examples
#' gradient_compute_pointwise(yy, Y, P)
#' @import Rfast
gradient_compute_pointwise = function(yy, Y, P){
  Y = c(t(rbind(yy,Y)))
  n = dim(P)[1]
  n2 = dim(P)[2]
  n_Y = length(Y)
  if ((n!=n2) | (2*n!=n_Y)){
    stop('Input dimensions do not match!')
  }

  YMat = matrix(Y, 2, n)
  YDiffDistSq = as.matrix(Dist(t(YMat))**2)
  YDiffDistSq_double = matrix(rep(YDiffDistSq[1,], each = 2), nrow = 2)
  P_double = matrix(rep(P[1,], each = 2), nrow = 2)
  deno = sum(1 / (1 + YDiffDistSq)) - n

  yy_double = matrix(rep(yy,n),2,n)
  yy_Y_diff = yy_double - YMat
  commonterm = yy_Y_diff / (1 + YDiffDistSq_double)
  I1 = 4 * P_double * commonterm
  I2 = (-1) * 4 * commonterm / deno / (1 + YDiffDistSq_double)

  G = rowsums(I1 + I2)

  return(G)
}

#' perturbation_score_compute_pointwise
#'
#' Given the perturbation direction and length for specified point, calculate the perturbation score.
#'
#' @param i Integer. Calculate the perturbation score for the i-th point.
#' @param X Matrix. Original Data.
#' @param Y Matrix. Embedding of original data.
#' @param perplexity Integer. Perplexity parameter to use, should be the same as when calculating the embedding Y.
#' @param dir_vec List. A list of perturbation directions.
#' @param length Numeric. Length of perturbation.
#' @param Ydist_sq Matrix. Distance matrix of embedding Y, can be calculated by as.matrix(dist(Y)**2).
#' @param initial_dims Integer. The number of dimensions that should be retained in the initial PCA step in t-SNE.
#' @param PCA_result PCA result of the original data. Can be calculated as prcomp(X).
#' @param P_unnormalized Matrix. Similarity matrix before normalization.
#' @param beta Vector. Vector of 1/(2sigma_i^2) in Gaussian kernel in input space in t-SNE.
#' @param approx Integer. 0: no approximation; 1: reuse the original PCA result; 2: with approximation to similarity matrix.
#' @return Perturbation score for the i-th point in the data.
#' @examples
#' perturbation_score_compute_pointwise(i, X, Y, perplexity, dir_vec, length, ...)
#' @export
#' @import RtsneWithP
perturbation_score_compute_pointwise = function(i, X, Y, perplexity, dir_vec, length, Ydist_sq = NULL, initial_dims = 50, PCA_result = NULL, P_unnormalized = NULL, beta = NULL, approx = 0){
  n = dim(X)[1]
  p = dim(X)[2]

  if (approx == 0){

    if (is.null(Ydist_sq)){
      Ydist_sq = as.matrix(dist(Y)**2)
    }

    YDistSqP1 = Ydist_sq
    YDistSqP1[2:n,2:n] = as.matrix(1 + Ydist_sq[-i,-i])
    YDistSqP1[1,1] = 1

    PScoreVec = numeric(6)
    dir_id = 1

    for (direction in dir_vec) {

      # Compute exact similarity matrix after perturbing X_i (P_new)
      X_new = rbind(X[i,] + direction * length, X[-i,])
      P_new = RtsneWithP::Rtsne(X_new,
                    initial_dims = initial_dims,
                    perplexity = perplexity,
                    theta = 0,
                    max_iter = 0,
                    Y_init = Y,
                    check_duplicates = FALSE)
      P_new = P_new$P

      best_value = Inf
      best_result = NULL

      # Take the embedding point before perturbation and the embedding points with the largest two similarity scores as initialization
      initializations_ind = order(P_new[1,], decreasing = TRUE)[1:2]
      initializations_ind = union(ifelse(initializations_ind == 1, i, ifelse(initializations_ind <= i, initializations_ind - 1, initializations_ind)), i)

      for (init_ind in initializations_ind) {

        # Find minimizer of LOO loss with Y[init_ind,] as initialization
        OPT = optim(Y[init_ind,],
                    fn = function(x){LOO_loss_fast(x, Y[-i,], YDistSqP1, P_new)},
                    gr = function(x){gradient_compute_pointwise(x, Y[-i,], P_new)},
                    method = 'BFGS', control = list(maxit = 100))

        if (OPT$value < best_value) {
          best_value = OPT$value
          best_result = OPT$par
        }
      }

      Y_new = best_result
      PScoreVec[dir_id] = sum((Y[i,] - Y_new) ^ 2) ^ 0.5
      dir_id = dir_id + 1
    }
    return(max(PScoreVec))

  }else if (approx == 1){

    if (is.null(Ydist_sq)){
      Ydist_sq = as.matrix(dist(Y)**2)
    }
    if (is.null(PCA_result)){
      stop('PCA result must be provided for if approx = 1!')
    }

    YDistSqP1 = Ydist_sq
    YDistSqP1[2:n,2:n] = as.matrix(1 + Ydist_sq[-i,-i])
    YDistSqP1[1,1] = 1

    PScoreVec = numeric(6)
    dir_id = 1

    for (direction in dir_vec) {

      # Compute exact similarity matrix after perturbing X_i (P_new)
      X_new = rbind(X[i,] + direction * length, X[-i,])
      P_new = RtsneWithP::Rtsne(X_new,
                    initial_dims = initial_dims,
                    perplexity = perplexity,
                    theta = 0,
                    max_iter = 0,
                    Y_init = Y,
                    check_duplicates = FALSE,
                    pca = FALSE)
      P_new = P_new$P

      best_value = Inf
      best_result = NULL

      # Take the embedding point before perturbation and the embedding points with the largest two similarity scores as initialization
      initializations_ind = order(P_new[1,], decreasing = TRUE)[1:2]
      initializations_ind = union(ifelse(initializations_ind == 1, i, ifelse(initializations_ind <= i, initializations_ind - 1, initializations_ind)), i)

      for (init_ind in initializations_ind) {

        # Find minimizer of LOO loss with Y[init_ind,] as initialization
        OPT = optim(Y[init_ind,],
                    fn = function(x){LOO_loss_fast(x, Y[-i,], YDistSqP1, P_new)},
                    gr = function(x){gradient_compute_pointwise(x, Y[-i,], P_new)},
                    method = 'BFGS', control = list(maxit = 100))

        if (OPT$value < best_value) {
          best_value = OPT$value
          best_result = OPT$par
        }
      }

      PScoreVec[dir_id] = sum((Y[i,] - best_result) ^ 2) ^ 0.5
      dir_id = dir_id + 1
    }
    return(max(PScoreVec))

  }else{

    if (is.null(Ydist_sq)){
      Ydist_sq = as.matrix(dist(Y)**2)
    }
    if (is.null(PCA_result)){
      stop('PCA result must be provided for if approx = 2!')
    }
    if (is.null(P_unnormalized)){
      stop('Unnormalized P must be provided for if approx = 2!')
    }
    if (is.null(beta)){
      stop('Beta must be provided for if approx = 2!')
    }

    beta = beta[-i]
    P_un = P_unnormalized
    P_unnormalized[2:n,2:n] = P_un[-i,-i]

    YDistSqP1 = Ydist_sq
    YDistSqP1[2:n,2:n] = as.matrix(1 + Ydist_sq[-i,-i])
    YDistSqP1[1,1] = 1

    PScoreVec = numeric(6)
    dir_id = 1

    for (direction in dir_vec) {

      # Compute approximated similarity matrix after perturbing X_i (P_new)
      X_new = rbind(X[i,] + direction * length, X[-i,])
      P_new = approxP(X[i,] + direction * length, X[-i,], beta, P_unnormalized, perplexity)

      # Take the embedding point before perturbation and the embedding points with the largest two similarity scores as initialization
      initializations_ind = order(P_new[1,], decreasing = TRUE)[1:2]
      initializations_ind = union(ifelse(initializations_ind == 1, i, ifelse(initializations_ind <= i, initializations_ind - 1, initializations_ind)), i)

      best_value = Inf
      best_result = NULL

      for (init_ind in initializations_ind) {

        # Find minimizer of LOO loss with Y[init_ind,] as initialization
        OPT = optim(Y[init_ind,],
                    fn = function(x){LOO_loss_fast(x, Y[-i,], YDistSqP1, P_new)},
                    gr = function(x){gradient_compute_pointwise(x, Y[-i,], P_new)},
                    method = 'BFGS', control = list(maxit = 100))

        if (OPT$value < best_value) {
          best_value = OPT$value
          best_result = OPT$par
        }
      }

      PScoreVec[dir_id] = sum((Y[i,] - best_result) ^ 2) ^ 0.5
      dir_id = dir_id + 1
    }
    return(max(PScoreVec))
  }
}

#' perturbation_score_compute
#'
#' Given the perturbation direction and length for specified point, calculate the perturbation scores.
#'
#' @param X Matrix. Original Data.
#' @param Y Matrix. Embedding of original data.
#' @param perplexity Integer. Perplexity parameter to use, should be the same as when calculating the embedding Y.
#' @param length Numeric. Length of perturbation.
#' @param approx Integer. 0: no approximation; 1: reuse the original PCA result; 2: with approximation to similarity matrix.
#' @param ind Vector. Specify indices of points for perturbation score calculation. (default: NULL (calculate for all points))
#' @param no.cores Integer. Number of cores for parallel computing. (default: NULL (no parallel computing))
#' @param initial_dims Integer. The number of dimensions that should be retained in the initial PCA step in t-SNE.
#' @param pca_center Logical. Should data be centered before pca is applied? (default: TRUE)
#' @param pca_scale Logical. Should data be scaled before pca is applied? (default: FALSE)
#' @return Perturbation score for the i-th point in the data.
#' @examples
#' perturbation_score_compute(X, Y, perplexity, length, approx)
#' @export
#' @import RtsneWithP
#' @import parallel
perturbation_score_compute = function(X, Y, perplexity, length, approx, ind = NULL, no.cores = NULL, initial_dims = 50, pca_center = TRUE, pca_scale = FALSE){
  if (approx){
    dir_length = min(initial_dims, dim(X)[2])
    PCA_x = prcomp(X, retx=TRUE, center = pca_center, scale. = pca_scale, rank. = dir_length)
    X = PCA_x$x
    PCA_x = prcomp(X)
  }else{
    dir_length = dim(X)[2]
    PCA_x = prcomp(X)
  }

  # calculate unnormalized similarity matrix and beta (1/sigma^2) for approximation
  P_unnormalized = NULL
  beta = NULL
  if (approx == 2){
    out = neMDBD::get_P_tsne(normalize_input(X), perplexity)
    P_unnormalized = out$P_unnormalized
    beta = out$beta
  }

  if (dim(X)[2] == 2){
    dir_vec = list(
      PCA_x$rotation[1:dir_length, 1],-PCA_x$rotation[1:dir_length, 1],
      PCA_x$rotation[1:dir_length, 2],-PCA_x$rotation[1:dir_length, 2]
    )
  }else{
    dir_vec = list(
      PCA_x$rotation[1:dir_length, 1],-PCA_x$rotation[1:dir_length, 1],
      PCA_x$rotation[1:dir_length, 2],-PCA_x$rotation[1:dir_length, 2],
      PCA_x$rotation[1:dir_length, 3],-PCA_x$rotation[1:dir_length, 3]
    )
  }

  # calculate distance matrix beforehand to save computation
  Ydist_sq = as.matrix(Dist(Y, square = TRUE))

  ind = if (is.null(ind)) c(1:dim(X)[1]) else ind

  if (is.null(no.cores)){
    pscore = sapply(ind, function(i){
      return(neMDBD::perturbation_score_compute_pointwise(i, X, Y, perplexity = perplexity, dir_vec = dir_vec, length = 1, Ydist_sq = Ydist_sq, initial_dims = 50, PCA_result = PCA_x, P_unnormalized = P_unnormalized, beta = beta, approx = approx))})
    return(pscore)
  }else{
    pscore = unlist(mclapply(ind,
                             function(i){return(neMDBD::perturbation_score_compute_pointwise(i, X, Y, perplexity = perplexity, dir_vec = dir_vec, length = 1, Ydist_sq = Ydist_sq, initial_dims = 50, PCA_result = PCA_x, P_unnormalized = P_unnormalized, beta = beta, approx = approx))},
                             mc.cores = no.cores))
    return(pscore)
  }
}

#####################################################################
#' get_P_tsne
#'
#' R function to calculate P matrix of t-SNE (https://github.com/jdonaldson/rtsne)
#'
#' @param X Pre-processed X maatrix
#' @param perplexity Perplexity of t-SNE algorithm.
#' @return P matrix, unnormalized P matrix, and beta (1/2sigma^2).
#' @examples
#' get_P_tsne(X, perplexity)
#' @import Rfast
get_P_tsne = function(X, perplexity) {
  .Hbeta = function(D, beta) {
    P = exp(-D * beta)
    sumP = sum(P)
    if (sumP == 0) {
      H = 0
      P = D * 0
    } else {
      H = log(sumP) + beta * sum(D %*% P) / sumP
      P = P / sumP
    }
    r = {
    }
    r$H = H
    r$P = P
    r
  }
  .Hbeta_unnormalized = function(D, beta) {
    P = exp(-D * beta)
    sumP = sum(P)
    if (sumP == 0) {
      H = 0
      P = D * 0
    } else {
      H = log(sumP) + beta * sum(D %*% P) / sumP
    }
    r = {
    }
    r$H = H
    r$P = P
    r
  }
  .x2p = function(X,
                  perplexity = 15,
                  tol = 1e-5) {
    if ('dist' %in% class(X)) {
      D = X
      n = attr(D, 'Size')
    } else{
      D = Dist(X) ** 2
      n = dim(X)[1]
    }

    D = as.matrix(D)
    P = matrix(0, n, n)
    P_unnormalized = matrix(0, n, n)
    beta = rep(1, n)
    logU = log(perplexity)

    for (i in 1:n) {
      betamin = -Inf
      betamax = Inf
      Di = D[i,-i]
      hbeta = .Hbeta(Di, beta[i])
      H = hbeta$H

      thisP = hbeta$P
      Hdiff = H - logU

      tries = 0


      while (abs(Hdiff) > tol && tries < 200) {
        if (Hdiff > 0) {
          betamin = beta[i]
          if (is.infinite(betamax))
            beta[i] = beta[i] * 2
          else
            beta[i] = (beta[i] + betamax) / 2
        } else{
          betamax = beta[i]
          if (is.infinite(betamin))
            beta[i] = beta[i] / 2
          else
            beta[i] = (beta[i] + betamin) / 2
        }

        hbeta = .Hbeta(Di, beta[i])
        H = hbeta$H
        thisP = hbeta$P
        Hdiff = H - logU
        tries = tries + 1
      }
      hbeta_unnormalized = .Hbeta_unnormalized(Di, beta[i])
      thisP_unnormalized = hbeta_unnormalized$P
      P[i, -i]  = thisP
      P_unnormalized[i, -i] = thisP_unnormalized
    }

    r = {
    }
    r$P = P
    r$P_unnormalized = P_unnormalized
    r$beta = beta
    sigma = sqrt(1 / beta)
    r
  }
  eps = 0#2^(-52)
  out = .x2p(X, perplexity, 1e-05)
  P = out$P
  beta = out$beta
  P = 0.5 * (P + t(P))
  P[P < eps] <- eps
  P = P / sum(P)
  return(list(P = P, P_unnormalized = out$P_unnormalized, beta = beta))
}

#' approxP
#'
#' R function to calculate approximated P matrix of t-SNE
#'
#' @param x_new The new x.
#' @param X The rest of X matrix.
#' @param beta beta for corresponding X.
#' @param P_unnormalized The unnormalized P matrix rearranged to put x_new as the first one.
#' @return approximated P matrix when changing x_new.
#' @examples
#' approxP(x_new, X, beta, P_unnormalized, perplexity)
#' @import Rfast
approxP = function(x_new, X, beta, P_unnormalized, perplexity){
  .Hbeta = function(D, beta) {
    P = exp(-D * beta)
    sumP = sum(P)
    if (sumP == 0) {
      H = 0
      P = D * 0
    } else {
      H = log(sumP) + beta * sum(D %*% P) / sumP
      P = P / sumP
    }
    r = {
    }
    r$H = H
    r$P = P
    r
  }
  .Hbeta_unnormalized = function(D, beta) {
    P = exp(-D * beta)
    sumP = sum(P)
    if (sumP == 0) {
      H = 0
      P = D * 0
    } else {
      H = log(sumP) + beta * sum(D %*% P) / sumP
    }
    r = {
    }
    r$H = H
    r$P = P
    r
  }


  X = normalize_input(rbind(x_new,X))
  n = dim(X)[1]
  p = dim(X)[2]
  Xdistsq = rowsums((X[-1,] - matrix(c(X[1,]), nrow = n-1, ncol = p, byrow = TRUE))**2)
  P_unnormalized[2:n,1] = exp(-Xdistsq * beta)

  logU = log(perplexity)
  betamin = -Inf
  betamax = Inf
  Di = Xdistsq
  beta1 = 1
  tol = 1e-5
  hbeta = .Hbeta(Di, beta1)
  H = hbeta$H

  thisP = hbeta$P
  Hdiff = H - logU

  tries = 0

  while (abs(Hdiff) > tol && tries < 200) {
    if (Hdiff > 0) {
      betamin = beta1
      if (is.infinite(betamax))
        beta1 = beta1 * 2
      else
        beta1 = (beta1 + betamax) / 2
    } else{
      betamax = beta1
      if (is.infinite(betamin))
        beta1 = beta1 / 2
      else
        beta1 = (beta1 + betamin) / 2
    }

    hbeta = .Hbeta(Di, beta1)
    H = hbeta$H
    thisP = hbeta$P
    Hdiff = H - logU
    tries = tries + 1
  }
  hbeta_unnormalized = .Hbeta_unnormalized(Di, beta1)
  thisP_unnormalized = hbeta_unnormalized$P
  P_unnormalized[1,2:n] = thisP_unnormalized
  P_new = P_unnormalized/rowsums(P_unnormalized)
  P_new = 0.5 * (P_new + t(P_new))
  P_new = P_new / sum(colsums(P_new))
  return(P_new)
}
