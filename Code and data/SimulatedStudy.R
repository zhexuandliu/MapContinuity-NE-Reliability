########################################
# Oct 31, 2024, Zhexuan Liu #############
# Simulations ##########################
########################################

library(ggplot2)
library(viridis)
library(neMDBD)
library(Rfast)
library(latex2exp)
library(clusterSim)

########################################
# 5GMM
########################################

# Original data, embedding, perturbation score (calculated by server) are stored in the following file.
load('./data/SimulatedStudy/5GMM_PScore_tsne_perplexity_60.RData')
px_label = ggplot(data = data.frame(X1 = X[,1], X2 = X[,2], label = label), aes(x = X1, y = X2, color = as.factor(label))) +
  geom_point(size = 2, show.legend = TRUE) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  ggtitle('Input Space with Label') + labs(color = 'Label')
py_label = ggplot(data = data.frame(X1 = Y[,1], X2 = Y[,2], label = label), aes(x = X1, y = X2, color = as.factor(label))) +
  geom_point(size = 1.5, show.legend = TRUE) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  ggtitle('Embedding with Label') + labs(color = 'Label') +  xlab('tSNE1') + ylab('tSNE2')
px_score = ggplot(data = data.frame(X1 = X[,1], X2 = X[,2], color = perturb_score),
                  aes(x = X1, y = X2, color = color)) +
  geom_point(size = 2) +
  scale_color_viridis() + labs(color = 'Perturbation\nScore') + ggtitle('Input Space with Perturbation Score')
py_score = ggplot(data = data.frame(X1 = Y[,1], X2 = Y[,2], color = perturb_score, point_size = ifelse(perturb_score > 8, 2.2, 1.5)),
                  aes(x = X1, y = X2, color = color, size = point_size)) +
  geom_point(show.legend = c(color = TRUE, size = FALSE)) +
  scale_size_identity() + scale_color_viridis() + labs(color = 'Perturbation\nScore') + ggtitle('Embedding with Perturbation Score') + xlab('tSNE1') + ylab('tSNE2')
ggsave(px_label, filename = './plots/SimulatedStudy/5GMM_pscore_x_label.png', width = 5, height = 4)
ggsave(py_label, filename = './plots/SimulatedStudy/5GMM_pscore_y_label.png', width = 5, height = 4)
ggsave(px_score, filename = './plots/SimulatedStudy/5GMM_pscore_x_score.png', width = 5, height = 4)
ggsave(py_score, filename = './plots/SimulatedStudy/5GMM_pscore_y_score.png', width = 5, height = 4)


#######################################
##### Entropy difference

# Input distance
mean_vec = list()
cov_vec = list()
cov_vec_diag = list()
mean_mat = matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
for (cl in unique(label)){
  mean_vec = c(mean_vec, list(colMeans(X[label == cl,])))
  mean_mat[label == cl,] = t(replicate(sum(label == cl),colMeans(X[label == cl,])))
  cov_vec = c(cov_vec, list(diag(apply(X[label == cl,], 2, var))))
}
COV = apply(X - mean_mat,2, var)
prob_X = sapply(1:length(mean_vec), function(i) dmvnorm(X[,!diag(cov_vec[[i]])<1e-10], mean_vec[[i]][!diag(cov_vec[[i]])<1e-10], cov_vec[[i]][!diag(cov_vec[[i]])<1e-10,!diag(cov_vec[[i]])<1e-10])) # prob
prob_prop_X = prob_X/rowsums(prob_X)
entropy_X = sapply(1:dim(X)[1], function(i) {-sum(prob_prop_X[i, prob_prop_X[i,]>1e-14] * log(prob_prop_X[i, prob_prop_X[i,]>1e-14]))})

# Output distance
mean_vec = list()
cov_vec = list()
mean_mat = matrix(0, nrow = dim(Y)[1], ncol = dim(Y)[2])
for (cl in unique(label)){
  mean_vec = c(mean_vec, list(colMeans(Y[label == cl,])))
  mean_mat[label == cl,] = t(replicate(sum(label == cl),colMeans(Y[label == cl,])))
  cov_vec = c(cov_vec, list(diag(apply(Y[label == cl,], 2, var))))
}
COV = apply(Y - mean_mat,2, var)
prob_Y = sapply(1:length(mean_vec), function(i) dmvnorm(Y, mean_vec[[i]], cov_vec[[i]])) # prob
prob_prop_Y = prob_Y/rowsums(prob_Y)
entropy_Y = sapply(1:dim(Y)[1], function(i) {-sum(prob_prop_Y[i, prob_prop_Y[i,]>1e-14] * log(prob_prop_Y[i, prob_prop_Y[i,]>1e-14]))})
entropy_difference = entropy_X - entropy_Y

px_entropy = ggplot(data = data.frame(X1 = X[,1], X2 = X[,2], color = entropy_difference), aes(x = X1, y = X2, color = color)) +
  geom_point(size = 2) +
  scale_color_viridis() + labs(color = 'Entropy\nDifference') + ggtitle('Input Space with Entropy Difference')
py_entropy = ggplot(data = data.frame(X1 = Y[,1], X2 = Y[,2], color = entropy_difference, point_size = ifelse(perturb_score > 8, 2.2, 1.5)), aes(x = X1, y = X2, color = color, size = point_size)) +
  geom_point(show.legend = c(color = TRUE, size = FALSE)) +
  scale_size_identity() + scale_color_viridis() + labs(color = 'Entropy\nDifference') + ggtitle('Embedding with Entropy Difference') + xlab('tSNE1') + ylab('tSNE2')
ggsave(px_entropy, filename = './plots/SimulatedStudy/5GMM_pscore_x_entropy.png', width = 5, height = 4)
ggsave(py_entropy, filename = './plots/SimulatedStudy/5GMM_pscore_y_entropy.png', width = 5, height = 4)


########################################
# 8GMM
########################################

perplexity_vec = c(5,50)
load(paste('./data/SimulatedStudy/8GMM_tsne_perplexity_',perplexity_vec[1],'.RData',sep=''))

sscore_list = vector("list", length = length(perplexity_vec))
plot_mat = data.frame(y1 = rep(0,length(perplexity_vec)*dim(X)[1]),
                      y2 = rep(0,length(perplexity_vec)*dim(X)[1]),
                      sscore = rep(0,length(perplexity_vec)*dim(X)[1]),
                      perplexity = rep(0,length(perplexity_vec)*dim(X)[1]),
                      label = rep(label,length(perplexity_vec)),
                      bi_label = rep(0,length(perplexity_vec)*dim(X)[1]))
pcol_list = vector('list', length = length(perplexity_vec))
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  load(paste('./data/SimulatedStudy/8GMM_tsne_perplexity_',perplexity,'.RData',sep=''))
  sscore_tmp = neMDBD::singularity_score_compute(Y,P)
  sscore_list[[i]] = sscore_tmp
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(1)] = Y[,1]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(2)] = Y[,2]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),3] = sscore_tmp
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),4] = perplexity
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),6] = (sscore_tmp > 2000)
}
options(scipen = -1)
my_labeller = function(labels) {return(paste("Perplexity", labels))}

p1 = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = label)) +
  geom_point(size = 0.3, show.legend = TRUE) +
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  xlab('tSNE1') + ylab('tSNE2') + labs(color = "Label") + guides(color = guide_legend(override.aes = list(size = 2)))

plot_mat$sscore_adj = plot_mat$sscore
plot_mat$sscore_adj[plot_mat$sscore>10000] = 10000
custom_labels = function(x) {ifelse(x >=10000, TeX('$\\geq$1e+04'), format(x))}
p2 = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = sscore_adj)) +
  geom_point(size = 0.3, show.legend = TRUE) +
  scale_color_viridis(direction = 1, trans = 'log10', name = "Singularity\nScore", labels = custom_labels) +
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  xlab('tSNE1') + ylab('tSNE2')

colors = c('#FF1F5B','#A0B1BA')
plot_mat$bi_label = factor(plot_mat$bi_label, levels = c('1', '0'))
p3 = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = as.factor(bi_label), size = ifelse(bi_label == 1, 1, 0.5))) +
  geom_point(show.legend = c(color = TRUE, size = FALSE)) +
  scale_size_identity()+
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  scale_color_manual(values = colors,name = "Dichotomized\nSingularity\nScore", labels = c('score\n> 2000', 'otherwise')) +
  guides(color = guide_legend(override.aes = list(size = 2))) + xlab('tSNE1') + ylab('tSNE2')

ggsave(p1,filename = './plots/SimulatedStudy/8GMM_plot1_compare2perps.png', width = 8, height = 10)
ggsave(p2,filename = './plots/SimulatedStudy/8GMM_plot1_compare2perps_sscore.png', width = 8, height = 10)
ggsave(p3,filename = './plots/SimulatedStudy/8GMM_plot1_compare2perps_bilabel.png', width = 8, height = 10)

#######################################
# test

# Perplexity 5
filename = paste('./data/SimulatedStudy/8GMM_tsne_perplexity_',5,'.RData',sep='')
load(filename)
kmeans_result <- kmeans(Y, centers = length(unique(label)), nstart = 100)
clusters <- kmeans_result$cluster

#### Davies-Bouldin Index
db_index <- index.DB(Y, clusters)$DB
print(paste("Davies-Bouldin Index:", db_index))
# [1] "Davies-Bouldin Index: 0.598224725466497"

#### within-cluster distance ratio
within_cluster_sumsq <- kmeans_result$tot.withinss
total_sumsq <- sum((Y - colMeans(Y))^2)
within_cluster_distance_ratio <- within_cluster_sumsq / total_sumsq
print(paste("Within-Cluster Distance Ratio:", within_cluster_distance_ratio))
# [1] "Within-Cluster Distance Ratio: 0.0480344447922012"

#### Wilks' lambda
fit = manova(Y ~ label, data = data.frame(Y=Y, label = label))
print(summary(fit, test = "Wilks"))
# Df     Wilks approx F num Df den Df
# label       7 0.0028423   2006.5     14   1582
# Residuals 792
# Pr(>F)
# label     < 2.2e-16 ***
#   Residuals
# ---
#   Signif. codes:
#   0 ‘***’ 1e-03 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Perplexity 50
filename = paste('./data/SimulatedStudy/8GMM_tsne_perplexity_',50,'.RData',sep='')
load(filename)
kmeans_result <- kmeans(Y, centers = length(unique(label)), nstart = 100)
clusters <- kmeans_result$cluster

#### Davies-Bouldin Index
db_index <- index.DB(Y, clusters)$DB
print(paste("Davies-Bouldin Index:", db_index))
# [1] "Davies-Bouldin Index: 0.303813334420334"

#### within-cluster distance ratio
within_cluster_sumsq <- kmeans_result$tot.withinss
total_sumsq <- sum((Y - colMeans(Y))^2)
within_cluster_distance_ratio <- within_cluster_sumsq / total_sumsq
print(paste("Within-Cluster Distance Ratio:", within_cluster_distance_ratio))
# [1] "Within-Cluster Distance Ratio: 0.00244015193998804"

#### Wilks' lambda
fit = manova(Y ~ label, data = data.frame(Y=Y, label = label))
print(summary(fit, test = "Wilks"))
# Df      Wilks approx F num Df den Df
# label       7 8.9514e-06    37656     14   1582
# Residuals 792
# Pr(>F)
# label     < 2.2e-16 ***
#   Residuals
# ---
#   Signif. codes:
#   0 ‘***’ 1e-03 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


########################################
# Swiss roll
########################################

perplexity_vec = c(5,40)
load(paste('./data/SimulatedStudy/swissroll_tsne_perplexity_',perplexity_vec[1],'.RData',sep=''))

sscore_list = vector("list", length = length(perplexity_vec))
plot_mat = data.frame(y1 = rep(0,length(perplexity_vec)*dim(X)[1]),
                      y2 = rep(0,length(perplexity_vec)*dim(X)[1]),
                      sscore = rep(0,length(perplexity_vec)*dim(X)[1]),
                      perplexity = rep(0,length(perplexity_vec)*dim(X)[1]),
                      label = rep(t,length(perplexity_vec)),
                      bi_label = rep(0,length(perplexity_vec)*dim(X)[1]))
pcol_list = vector('list', length = length(perplexity_vec))
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  load(paste('./data/SimulatedStudy/swissroll_tsne_perplexity_',perplexity,'.RData',sep=''))
  sscore_tmp = neMDBD::singularity_score_compute(Y,P)
  sscore_list[[i]] = sscore_tmp
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(1)] = Y[,1]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(2)] = Y[,2]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),3] = sscore_tmp
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),4] = perplexity
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),6] = (sscore_tmp > 4000)
}
options(scipen = -1)
my_labeller = function(labels) {return(paste("Perplexity", labels))}

p1 = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = label)) +
  geom_point(size = 0.3, show.legend = TRUE) +
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  xlab('tSNE1') + ylab('tSNE2') + labs(color = 'Spiral\nAngle')

plot_mat$sscore_adj = plot_mat$sscore
plot_mat$sscore_adj[plot_mat$sscore>10000] = 10000
custom_labels = function(x) {ifelse(x >=10000, TeX('$\\geq$1e+04'), format(x))}
p2 = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = sscore_adj)) +
  geom_point(size = 0.3, show.legend = TRUE) +
  scale_color_viridis(direction = 1, trans = 'log10', name = "Singularity\nScore", labels = custom_labels) +
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  xlab('tSNE1') + ylab('tSNE2')

colors = c('#FF1F5B','#A0B1BA')
plot_mat$bi_label = factor(plot_mat$bi_label, levels = c('1', '0'))
p3 = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = as.factor(bi_label), size = ifelse(bi_label == 1, 1, 0.5))) +
  geom_point(show.legend = c(color = TRUE, size = FALSE)) +
  scale_size_identity()+
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  scale_color_manual(values = colors,name = "Dichotomized\nSingularity\nScore", labels = c('score\n> 4000', 'otherwise')) +
  guides(color = guide_legend(override.aes = list(size = 2))) + xlab('tSNE1') + ylab('tSNE2')

ggsave(p1,filename = './plots/SimulatedStudy/Swissroll_plot1_compare2perps.png', width = 8, height = 10)
ggsave(p2,filename = './plots/SimulatedStudy/Swissroll_plot1_compare2perps_sscore.png', width = 8, height = 10)
ggsave(p3,filename = './plots/SimulatedStudy/Swissroll_plot1_compare2perps_bilabel.png', width = 8, height = 10)
