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
# OI and FI discontinuity in UMAP
library(uwot)
library(ggplot2)
library(MGMM)

set.seed(4)
X = MGMM::rGMM(500, d = 2, k = 2, means = list(c(2, 0), c(-2, 0)), covs = diag(2))
label = factor(rownames(X))
colors <- c('#FF1F5B','#009ADE')
p_x = ggplot() +  
  geom_point(data = data.frame(x = X[,1], y = X[,2], label = label), 
             aes(x=x,y=y, color=label), size = 1.2) + 
  xlab('X1') + 
  ylab('X2') + 
  theme(legend.position = 'None') +
  ggtitle('Input Space with Labels')+
  scale_color_manual(values = colors)+ 
  theme_minimal()+
  scale_color_manual(values = colors,name = "Label")+
  guides(color = guide_legend(override.aes = list(size = 2)))

Y = umap(X, n_neighbors = 5)
p_y_5 = ggplot() +  
  geom_point(data = data.frame(x = Y[,1], y = Y[,2], label = label), 
             aes(x=x,y=y, color=label), size = 1.2) + 
  xlab('UMAP1') + 
  ylab('UMAP2') + 
  theme(legend.position = 'None') +
  ggtitle('Number of neighbors 5') + 
  theme_minimal()+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16),
    axis.line = element_blank(),
    legend.position = "none"
  ) + 
  scale_color_manual(values = colors,name = "Label")

Y = umap(X, n_neighbors = 50)
p_y_50 = ggplot() +  
  geom_point(data = data.frame(x = Y[,1], y = Y[,2], label = label), 
             aes(x=x,y=y, color=label), size = 1.2) + 
  xlab('UMAP1') + 
  ylab('UMAP2') + 
  theme(legend.position = 'None') +
  ggtitle('Number of neighbors 50') + 
  theme_minimal()+
  scale_color_manual(values = colors,name = "Label")

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

#######################################
# Supplementary Fig. 2
## a
library(ggplot2)
source('./data/2GMM/disconuity plot.R')
load('./data/2GMM/2d_2GMM_perplexity5.RData')
perplexity = 5
mean_1 = colmeans(X[label == 1,])
mean_2 = colmeans(X[label == 2,])
P1 = plot_interpolate(X, perplexity, Y, mean_2, mean_1, label, 999, c(500), 1,size_all = 1.2, size_min = 7, method = 'tsne', dens = 100, thr = 100) + ggtitle(paste('Perplexity',perplexity,sep=' '))
load('./data/2GMM/2d_2GMM_perplexity20.RData')
perplexity = 20
mean_1 = colmeans(X[label == 1,])
mean_2 = colmeans(X[label == 2,])
P2 = plot_interpolate(X, perplexity, Y, mean_2, mean_1, label, 999, c(500), 1,size_all = 1.2, size_min = 7, method = 'tsne', dens = 100, thr = 100)+ ggtitle(paste('Perplexity',perplexity,sep=' '))
load('./data/2GMM/2d_2GMM_perplexity40.RData')
perplexity = 40
mean_1 = colmeans(X[label == 1,])
mean_2 = colmeans(X[label == 2,])
P3 = plot_interpolate(X, perplexity, Y, mean_2, mean_1, label, 999, c(500), 1,size_all = 1.2, size_min = 7, method = 'tsne', dens = 100, thr = 100)+ ggtitle(paste('Perplexity',perplexity,sep=' '))
load('./data/2GMM/2d_2GMM_perplexity60.RData')
perplexity = 60
mean_1 = colmeans(X[label == 1,])
mean_2 = colmeans(X[label == 2,])
P4 = plot_interpolate(X, perplexity, Y, mean_2, mean_1, label, 999, c(500), 1,size_all = 1.2, size_min = 7, method = 'tsne', dens = 100, thr = 100)+ ggtitle(paste('Perplexity',perplexity,sep=' '))
load('./data/2GMM/2d_2GMM_perplexity80.RData')
perplexity = 80
mean_1 = colmeans(X[label == 1,])
mean_2 = colmeans(X[label == 2,])
P5 = plot_interpolate(X, perplexity, Y, mean_2, mean_1, label, 999, c(500), 1,size_all = 1.2, size_min = 7, method = 'tsne', dens = 100, thr = 100)+ ggtitle(paste('Perplexity',perplexity,sep=' '))
load('./data/2GMM/2d_2GMM_perplexity100.RData')
perplexity = 100
mean_1 = colmeans(X[label == 1,])
mean_2 = colmeans(X[label == 2,])
P6 = plot_interpolate(X, perplexity, Y, mean_2, mean_1, label, 999, c(500), 1,size_all = 1.2, size_min = 7, method = 'tsne', dens = 100, thr = 100)+ ggtitle(paste('Perplexity',perplexity,sep=' '))

## b
plot_list = vector("list", 6)
perplexity_vec = c(5,20,40,60,80,100)
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  filename1 = paste('./data/2GMM/2d_2GMM_perplexity',perplexity,'.RData', sep = '')
  filename2 = paste('./data/2GMM/coord_min_perplexity',perplexity,'.RData', sep = '')
  load(filename1)
  load(filename2)
  plot_list[[i]] = ggplot(data = data.frame(x = Y[,1], y = Y[,2], label = label), aes(x=x,y=y)) + 
    geom_point(aes(color = label), size = 0.5, alpha = 0.4)+ geom_path(data = data.frame(x = coord_min[,1], y = coord_min[,2]),aes(x= x, y = y),linetype = 'dotted', size = 0.6) + geom_point(data = data.frame(x = coord_min[,1], y = coord_min[,2]),aes(x= x, y = y), color = 'orange', size = 3, shape = 17)+
    scale_color_manual(values = c('#FF1F5B', '#009ADE')) + 
    xlab('tSNE1') + ylab('tSNE2') +
    theme_minimal()+ ggtitle(paste('Perplexity',perplexity,sep=' '))
}

## c
perplexity_vec = c(5,20,40,60,80,100)
load('./data/2GMM/2d_2GMM_perplexity5.RData')
eigen_score_list = vector("list", length = length(perplexity_vec))
plot_mat = data.frame(y1 = rep(0,length(perplexity_vec)*dim(X)[1]),
                      y2 = rep(0,length(perplexity_vec)*dim(X)[1]),
                      score_eigen = rep(0,length(perplexity_vec)*dim(X)[1]),
                      perplexity = rep(0,length(perplexity_vec)*dim(X)[1]),
                      label = rep(0,dim(X)[1]*length(perplexity_vec)),
                      bi_label = rep(0,length(perplexity_vec)*dim(X)[1]),
                      theta = rep(0,length(perplexity_vec)*dim(X)[1]))
pcol_list = vector('list', length = length(perplexity_vec))
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  filename = paste('./data/2GMM/2d_2GMM_perplexity',perplexity,'.RData',sep='')
  load(filename)
  singularity_score = neMDBD::singularity_score_compute(Y, P)
  eigen_score_list[[i]] = singularity_score
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(1)] = Y[,1]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(2)] = Y[,2]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),3] = singularity_score
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),4] = perplexity
  # plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),6] = (1/singularity_score > quantile(1/singularity_score, 0.95)) # binarize by quantile
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),6] = (singularity_score > 3000) # binarize by threshold
}
options(scipen = -1)
my_labeller <- function(labels) {
  return(paste("Perplexity", labels))
}
colors <- c('#FF1F5B','#A0B1BA')
plot_mat$bi_label <- factor(plot_mat$bi_label, levels = c('1', '0'))
p3 <- ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = as.factor(bi_label))) +
  geom_point(size = 0.8, show.legend = TRUE)  +
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  xlab('tSNE1') + ylab('tSNE2') +
  theme(strip.text = element_text(colour = "black", face = "bold"),
        strip.background.x = element_rect(fill = "white"),
        strip.background.y = element_rect(fill = "white"))+
  # scale_color_manual(values = colors,name = "Label", labels = c('> 95%\nquantile', 'otherwise'))
  scale_color_manual(values = colors,name = "Dichotomized\nSingularity\nScore", labels = c('score\n> 3000', 'otherwise')) + 
  theme_minimal()+ 
  guides(color = guide_legend(override.aes = list(size = 2)))

#######################################
# Supplementary Fig. 9
load("./data/2GMM/2d_2GMM_perplexity50.RData")
load("./data/2GMM/perturbation_score_perplexity50.RData")
variance_scores = c(read.table("./data/2GMM/variance_scores_perplexity50.txt", header = FALSE, sep = "\t"))$V1
df <- data.frame(x = Y[,1], y = Y[,2], color_value = variance_scores)

# Plot
p1 = ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color = color_value), size = 1.2) +
  scale_color_gradientn(colors = c("blue", "lightblue", "white", "lightpink", "red"),
                        values = scales::rescale(c(0.01, 0.05, 0.08, 0.1, 0.12)/0.12*0.132),
                        limits = c(0, 0.132))+
  theme_minimal() +
  theme(legend.key.height = unit(1.5, "cm")) +
  labs(color = "variance\nscore") + 
  ggtitle('t-SNE Embedding with scores by DynamicViz')+ xlab('tSNE1') + ylab('tSNE2')
p2 = ggplot() +
  geom_point(data = data.frame(
    x = Y[, 1],
    y = Y[, 2],
    score = perturbation_score
  ),
  aes(x = x, y = y, color = score)) +
  viridis::scale_color_viridis(direction = 1,
                               name = "Perturbation\nScore") + 
  ggtitle('t-SNE Embedding with Perturbation Score')+
  xlab('tSNE1') + ylab('tSNE2')+
  theme_minimal()