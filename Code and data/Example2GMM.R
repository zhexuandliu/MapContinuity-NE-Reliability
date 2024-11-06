########################################
# Oct 31, 2024, Zhexuan Liu #############
# Two components Gaussian mixture data #
########################################

# load packages
library(neMDBD)
library(RtsneWithP)
library(ggplot2)
library(viridis)
library(Rfast)
library(parallel)

# generate data
set.seed(4)
X = MGMM::rGMM(500, d = 2, k = 2, means = list(c(2, 0), c(-2, 0)), covs = diag(2))
label = factor(rownames(X))


#######################################
# Observe OI and FI discontinuity
#######################################

# OI discontinuity when perplexity is large
perplexity = 50
tsne.out = Rtsne(X, perplexity = perplexity, theta = 0, max_iter = 1000, Y_init = X)
P = tsne.out$P
Y = tsne.out$Y
p_x = ggplot() +
  geom_point(data = data.frame(X1 = X[, 1], X2 = X[, 2], label = label), aes(x = X1, y = X2, color = label)) +
  scale_color_manual(values = c('#FF1F5B','#009ADE'))
p_y = ggplot() +
  geom_point(data = data.frame(tSNE1 = Y[, 1], tSNE2 = Y[, 2], label = label), aes(x = tSNE1, y = tSNE2, color = label)) +
  scale_color_manual(values = c('#FF1F5B','#009ADE'))
ggsave(p_x, file = './plots/2GMM/2GMM_original.png', width = 5, height = 4)
ggsave(p_y, file = './plots/2GMM/2GMM_embedding_perplexity50.png', width = 5, height = 4)

# FI discontinuity when perplexity is small
perplexity = 5
tsne.out = Rtsne(X, perplexity = perplexity, theta = 0, max_iter = 1000, Y_init = X)
P = tsne.out$P
Y = tsne.out$Y
p_x = ggplot() +
  geom_point(data = data.frame(X1 = X[, 1], X2 = X[, 2], label = label), aes(x = X1, y = X2, color = label)) +
  scale_color_manual(values = c('#FF1F5B','#009ADE'))
p_y = ggplot() +
  geom_point(data = data.frame(tSNE1 = Y[, 1], tSNE2 = Y[, 2], label = label), aes(x = tSNE1, y = tSNE2, color = label)) +
  scale_color_manual(values = c('#FF1F5B','#009ADE'))
ggsave(p_y, file = './plots/2GMM/2GMM_embedding_perplexity5.png', width = 5, height = 4)


#######################################
# Contour plot of LOO loss
#######################################

# function to plot:
plot_contour_LOO_loss = function(X, Y, perplexity, label, x_new, thr = Inf) {
  ## Plot the contour of LOO loss w.r.t. x_new.
  # X           Matrix. Original data matrix.
  # Y           Matrix. The embedding of X.
  # perplexity  Interger. The perplexity used when calculating Y.
  # label       Vector. The labels of X.
  # x_new       Vector. The new point in the original space.
  # thr         Numeric. Threshold the LOO loss when plotting to focus on the hyperbolic structure.
  
  # create plot data frame, store the LOO loss value at embedding point (x,y) in z
  xr = c(min(Y[,1]), max(Y[,1])) # margin for x axis
  yr = c(min(Y[,2]), max(Y[,2])) # margin for y axis
  dens = 100 # grid density
  xgrid = c(0:dens) / dens * (xr[2] - xr[1]) + xr[1]
  ygrid = c(0:dens) / dens * (yr[2] - yr[1]) + yr[1]
  plotdf = data.frame(x = c(rep(xgrid, each = length(ygrid))), 
                      y = c(rep(ygrid, length(xgrid))), 
                      z = c(rep(0, length(xgrid) * length(ygrid))))
  
  X_new = rbind(x_new, X)
  labelnew = c('new', label)
  P_new = RtsneWithP::Rtsne(X_new, perplexity = perplexity, theta = 0, max_iter = 0, Y_init = rbind(c(0,0),Y), check_duplicates = FALSE)$P
  
  YDistSqP1 = as.matrix(Dist(rbind(c(0,0), Y), square = TRUE)) + 1
  plotdf$z = unlist(mclapply(1:dim(plotdf)[1],
                             function(i) return(LOO_loss_fast(c(plotdf$x[i], plotdf$y[i]), Y, YDistSqP1, P_new)),
                             mc.cores = detectCores() - 1))
  
  # plot
  colnames(Y) = c('x','y')
  p_contour = ggplot(data = plotdf) +
    geom_contour(aes(x = x, y = y, z = ifelse(z>thr, thr, z), colour = after_stat(level)), bins = 8, size = 1) +
    geom_point(data = data.frame(Y[label == unique(label)[1],],label[label == unique(label)[1]]), aes(x = x, y = y), color = '#009ADE', size = 1.2, alpha = 0.25)+
    geom_point(data = data.frame(Y[label == unique(label)[2],],label[label == unique(label)[2]]), aes(x = x, y = y), color = '#FF1F5B', size = 1.2, alpha = 0.25)+
    geom_point(data = plotdf[which.min(plotdf$z), ], aes(x = x, y = y), color = 'orange', size = 7, shape = 17) + 
    scale_color_viridis_c(direction = 1) + xlab('tSNE1') + ylab('tSNE2') + guides(color = guide_colorbar(title = "LOO Loss"))
  
  return(p_contour)
}

# Perplexity 50
perplexity = 50
load('./data/2GMM/2d_2GMM_perplexity50.RData')
mean_1 = colmeans(X[label == 1,])
mean_2 = colmeans(X[label == 2,])

x_new = mean_2
P1 = plot_contour_LOO_loss(X, Y, perplexity, label, x_new, thr = 8.462)
x_new = 0.47 * mean_1 + 0.53 * mean_2
P2 = plot_contour_LOO_loss(X, Y, perplexity, label, x_new, thr = 8.462)
x_new = 0.48 * mean_1 + 0.52 * mean_2
P3 = plot_contour_LOO_loss(X, Y, perplexity, label, x_new, thr = 8.462)
x_new = mean_1
P4 = plot_contour_LOO_loss(X, Y, perplexity, label, x_new, thr = 8.462)

limit = c(8.4552, 8.462)
P = egg::ggarrange(P1 + theme(legend.position = 'none') + scale_color_viridis_c(limits = limit),
                   P2 + theme(legend.position = 'none') + scale_color_viridis_c(limits = limit),
                   P3 + theme(legend.position = 'none') + scale_color_viridis_c(limits = limit),
                   P4 + scale_color_viridis_c(limits = limit), nrow = 1)
ggsave(P, file = './plots/2GMM/2GMM_contour_perplexity50.png', width = 20, height = 4)


# Perplexity 5
perplexity = 5
load('./data/2GMM/2d_2GMM_perplexity5.RData')
mean_1 = colmeans(X[label == 1,])
mean_2 = colmeans(X[label == 2,])

x_new = 0.2 * mean_1 + 0.8 * mean_2
P1 = plot_contour_LOO_loss(X, Y, perplexity, label, x_new)
x_new = 0.4 * mean_1 + 0.6 * mean_2
P2 = plot_contour_LOO_loss(X, Y, perplexity, label, x_new)
x_new = 0.6 * mean_1 + 0.4 * mean_2
P3 = plot_contour_LOO_loss(X, Y, perplexity, label, x_new)
x_new = 0.8 * mean_1 + 0.2 * mean_2
P4 = plot_contour_LOO_loss(X, Y, perplexity, label, x_new)

limit = c(5.8605, 5.9090)
P = egg::ggarrange(P1 + theme(legend.position = 'none') + scale_color_viridis_c(limits = limit),
                   P2 + theme(legend.position = 'none') + scale_color_viridis_c(limits = limit),
                   P3 + theme(legend.position = 'none') + scale_color_viridis_c(limits = limit),
                   P4 + scale_color_viridis_c(limits = limit), nrow = 1)
ggsave(P, file = './plots/2GMM/2GMM_contour_perplexity5.png', width = 20, height = 4)
