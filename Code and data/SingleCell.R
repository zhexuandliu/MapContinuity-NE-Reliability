########################################
# Oct 31, 2024, Zhexuan Liu #############
# Single-cell data #####################
########################################

library(ggplot2)
library(viridis)
library(neMDBD)

########################################
# Mouse brain embryo differentiation
########################################

#######################################
##### compare small perplexity and large perplexity

perplexity_vec = c(4,25)
load('./data/SingleCell/hayashi.RData')
load(paste('./data/SingleCell/hayashi_tsne_perplexity_',perplexity_vec[1],'.RData',sep=''))
sscore_list = vector("list", length = length(perplexity_vec))
plot_mat = data.frame(y1 = rep(0,length(perplexity_vec)*dim(X)[1]),
                      y2 = rep(0,length(perplexity_vec)*dim(X)[1]),
                      sscore = rep(0,length(perplexity_vec)*dim(X)[1]),
                      perplexity = rep(0,length(perplexity_vec)*dim(X)[1]),
                      label = rep(as.factor(ct.X1),length(perplexity_vec)),
                      bi_label = rep(0,length(perplexity_vec)*dim(X)[1]))
pcol_list = vector('list', length = length(perplexity_vec))
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  filename = paste('./data/SingleCell/hayashi_tsne_perplexity_',perplexity,'.RData',sep='')
  load(filename)
  sscore = neMDBD::singularity_score_compute(Y,P)
  sscore_list[[i]] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(1)] = Y[,1]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(2)] = Y[,2]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),3] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),4] = perplexity
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),6] = (sscore > 1800)
}
options(scipen = -1)
my_labeller = function(labels) {return(paste("Perplexity", labels))}

# Create the plot of singularity score (dichotomized)
colors = c('#FF1F5B','#A0B1BA')
plot_mat$bi_label = factor(plot_mat$bi_label, levels = c('1', '0'))
p_dicho = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = as.factor(bi_label), size = as.factor(bi_label))) +
  geom_point(show.legend = TRUE)+
  scale_size_manual(values = c("1" = 1, "0" = 1))+  
  guides(size = "none")  +
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  xlab('tSNE1') + ylab('tSNE2') +
  scale_color_manual(values = colors,name = "Dichotomized\nSingularity\nScore", labels = c('score\n> 1800', 'otherwise')) + 
  guides(color = guide_legend(override.aes = list(size = 2)))

# plot by label
label_time = as.numeric(plot_mat$label)
label_time[label == 1] = '00 h'
label_time[label == 2] = '12 h'
label_time[label == 3] = '24 h'
label_time[label == 4] = '48 h'
label_time[label == 5] = '72 h'
plot_mat$label_time= as.factor(label_time)
p_label = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = label_time)) +
  geom_point(size = 1.2, show.legend = TRUE) +
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  xlab('tSNE1') + ylab('tSNE2') +
  labs(color = "Time\nPoint") + 
  guides(color = guide_legend(override.aes = list(size = 2)))

P = egg::ggarrange(p_label, p_dicho, nrow = 1)
ggsave(P, filename = './plots/SingleCell/embryo/hayashi_plot.png', width = 10, height = 8)

# degree of FI discontinuity decreases in perplexity

perplexity_vec = c(seq(2,40,1))
load('./data/SingleCell/hayashi_tsne_perplexity_2.RData')
sscore_list = vector("list", length = length(perplexity_vec))
plot_mat = data.frame(y1 = rep(0,length(perplexity_vec)*dim(X)[1]),
                      y2 = rep(0,length(perplexity_vec)*dim(X)[1]),
                      sscore = rep(0,length(perplexity_vec)*dim(X)[1]),
                      perplexity = rep(0,length(perplexity_vec)*dim(X)[1]))
pcol_list = vector('list', length = length(perplexity_vec))
# the following for-loop may take a few minutes to execute
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  load(paste('./data/SingleCell/hayashi_tsne_perplexity_',perplexity,'.RData', sep = ''))
  sscore = neMDBD::singularity_score_compute(Y,P)#
  sscore_list[[i]] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(1,2)] = Y
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),3] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),4] = perplexity
}

q = 0.05
mean_1_sscore_q = numeric(length(perplexity_vec))
for (i in 1:length(perplexity_vec)){
  sscore = sscore_list[[i]]
  mean_1_sscore_q[[i]] = mean(sscore[sscore>quantile(sscore, 1-q)])
}
p_sscore = ggplot(data = data.frame(mean = mean_1_sscore_q, perplexity = perplexity_vec)) +
  geom_point(aes(x = perplexity, y = mean)) + ggtitle('Degree of FI Discontinuity')+ theme_classic()+
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.5), 
    panel.grid.minor = element_line(color = "grey95", size = 0.25) 
  ) + ylab('Mean of top 5% singularity score') + xlab('Perplexity') + 
  geom_point(data = data.frame(x = c(4,25), y = c(mean_1_sscore_q[c(3,24)])), aes(x = x, y = y), size = 2, colour='#FF1F5B')

#######################################
##### neighborhood preservation increases in perplexity

get_cor_score = function(X, Y, dist_X, nn_mat) {
  dist_Y = as.matrix(dist(Y))
  score_vec = sapply(c(1:dim(X)[1]), function(i) {
    nn_ind = nn_mat[i,]
    return(cor(c(dist_X[i,nn_ind]),c(dist_Y[i,nn_ind])))
  })
  return(score_vec)
}

perplexity_vec = c(seq(2,40,1))
load('./data/SingleCell/hayashi_tsne_perplexity_5.RData')
cor_score_vec = c(rep(0, length(perplexity_vec)))
pca_result = prcomp(X, rank = 50)
dist_X = as.matrix(dist(pca_result$x[,1:50]))
knn_result = FNN::get.knn(pca_result$x[,1:50], round(dim(X)[1]/5))
nn_mat = knn_result$nn.index
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  filename = paste('./data/SingleCell/hayashi_tsne_perplexity_',perplexity,'.RData', sep = '')
  load(filename)
  cor_score_vec[i] = median(get_cor_score(pca_result$x[,1:50], Y, dist_X, nn_mat))
}

p_nnpre = ggplot(data = data.frame(score = cor_score_vec, perplexity = perplexity_vec)) +
  geom_point(aes(x = perplexity, y = score)) + ggtitle('Neighborhood Preservation')+ theme_classic()+
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.5), 
    panel.grid.minor = element_line(color = "grey95", size = 0.25) 
  ) + ylab('Correlation of NN distances') + xlab('Perplexity') + 
  geom_point(data = data.frame(x = c(4,25), y = c(cor_score_vec[c(3,24)])), aes(x = x, y = y), size = 2, colour='#FF1F5B')

P = egg::ggarrange(p_sscore, p_nnpre, ncol = 1)
ggsave(P, filename = './plots/SingleCell/embryo/hayashi_choose_perplexity.png', width = 5, height = 8)


#######################################
##### tests

## Spearson's test
load('./data/SingleCell/hayashi_tsne_perplexity_4.RData')
singularity_score = neMDBD::singularity_score_compute(Y, P)
# output mean list for each group
mean_vec = list()
mean_mat = matrix(0, nrow = dim(Y)[1], ncol = dim(Y)[2])
label_unique = unique(label)
for (cl in unique(label)){
  mean_vec = c(mean_vec, list(colMeans(Y[label == cl,])))
  mean_mat[label == cl,] = t(replicate(sum(label == cl),colMeans(Y[label == cl,])))
}
# output distance to respective group center
dist2center = sqrt(rowSums((Y-mean_mat)**2))
# perform Spearman rank test in R
p_values = numeric()
for (cl in label_unique){
  test_result = cor.test(singularity_score[label == cl], dist2center[label == cl], method = "spearman")
  p_values = c(p_values, test_result$p.value)
}
print('Perplexity 4: p values are (may consider Holm-Bonferroni correction)')
print(p_values)
print(p.adjust(p_values, method = "holm"))

load('./data/SingleCell/hayashi_tsne_perplexity_25.RData')
singularity_score = neMDBD::singularity_score_compute(Y, P)
# output mean list for each group
mean_vec = list()
mean_mat = matrix(0, nrow = dim(Y)[1], ncol = dim(Y)[2])
label_unique = unique(label)
for (cl in unique(label)){
  mean_vec = c(mean_vec, list(colMeans(Y[label == cl,])))
  mean_mat[label == cl,] = t(replicate(sum(label == cl),colMeans(Y[label == cl,])))
}
# output distance to respective group center
dist2center = sqrt(rowSums((Y-mean_mat)**2))
# perform Spearman rank test in R
p_values = numeric()
for (cl in label_unique){
  test_result = cor.test(singularity_score[label == cl], dist2center[label == cl], method = "spearman")
  p_values = c(p_values, test_result$p.value)
}
print('Perplexity 25: p values are (may consider Holm-Bonferroni correction)')
print(p_values)
print(p.adjust(p_values, method = "holm"))

## F test for local regression
load('./data/SingleCell/hayashi_tsne_perplexity_4.RData')
singularity_score = neMDBD::singularity_score_compute(Y, P)
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = singularity_score, one = rep(1, dim(Y)[1]))
fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
fit_reduced = loess(score ~ -y1-y2, span = 1, degree = 0, data = data)
print(anova(fit_full,fit_reduced))
# Model 1: loess(formula = score ~ y1 * y2, data = data, span = 0.5, degree = 1)
# Model 2: loess(formula = score ~ -y1 - y2, data = data, span = 1, degree = 0)
# 
# Analysis of Variance:   denominator df 412.12
# 
# ENP        RSS F-value Pr(>F)
# [1,] 6.54 3.1085e+10               
# [2,] 1.94 3.1775e+10  1.1036 0.3598

load('./data/SingleCell/hayashi_tsne_perplexity_25.RData')
singularity_score = neMDBD::singularity_score_compute(Y, P)
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = singularity_score, one = rep(1, dim(Y)[1]))
fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
fit_reduced = loess(score ~ -y1-y2, span = 1, degree = 0, data = data)
print(anova(fit_full,fit_reduced))
# Model 1: loess(formula = score ~ y1 * y2, data = data, span = 0.5, degree = 1)
# Model 2: loess(formula = score ~ -y1 - y2, data = data, span = 1, degree = 0)
# 
# Analysis of Variance:   denominator df 411.11
# 
# ENP        RSS F-value    Pr(>F)    
# [1,] 7.19 1.3868e+09                      
# [2,] 1.61 1.4946e+09   3.272 9.709e-04 ***
#   ---
#   Signif. codes:  
#   0 ‘***’ 1e-03 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Permutation test for local regression
load('./data/SingleCell/hayashi_tsne_perplexity_4.RData')
singularity_score = neMDBD::singularity_score_compute(Y, P)
set.seed(1)
n = 10000
F_vec = numeric(n)
RSS_reduced = (dim(Y)[1]-1) * var(singularity_score)
for (i in c(1:n)){
  data = data.frame(y1 = Y[,1], y2 = Y[,2], score = sample(singularity_score))
  fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
  RSS_full = sum(fit_full$residuals**2)
  F_vec[i] = (RSS_reduced - RSS_full) / RSS_full
}
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = (singularity_score))
fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
RSS_full = sum(fit_full$residuals**2)
F_stat = (RSS_reduced - RSS_full) / RSS_full
print(sum(F_stat<F_vec)/length(F_vec))
# [1] 0.3195

load('./data/SingleCell/hayashi_tsne_perplexity_25.RData')
singularity_score = neMDBD::singularity_score_compute(Y, P)
set.seed(1)
n = 10000
F_vec = numeric(n)
RSS_reduced = (dim(Y)[1]-1) * var(singularity_score)
for (i in c(1:n)){
  data = data.frame(y1 = Y[,1], y2 = Y[,2], score = sample(singularity_score))
  fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
  RSS_full = sum(fit_full$residuals**2)
  F_vec[i] = (RSS_reduced - RSS_full) / RSS_full
}
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = (singularity_score))
fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
RSS_full = sum(fit_full$residuals**2)
F_stat = (RSS_reduced - RSS_full) / RSS_full
print(sum(F_stat<F_vec)/length(F_vec))
# [1] 0.0074

########################################
# 10X geno brain
########################################

#######################################
##### compare small perplexity and large perplexity

perplexity_vec = c(5,95)
load('./data/SingleCell/brain.RData')
X = t(data.X1)
load(paste('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_',perplexity_vec[1],'.RData',sep=''))
sscore_list = vector("list", length = length(perplexity_vec))
plot_mat = data.frame(y1 = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      y2 = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      sscore = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      perplexity = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      label = rep(as.factor(ct.X1),length(perplexity_vec)),
                      bi_label = rep(0,length(perplexity_vec)*dim(Y)[1]))
pcol_list = vector('list', length = length(perplexity_vec))
# the following for-loop may take a few minutes to execute
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  filename = paste('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_',perplexity,'_Y.RData',sep='')
  load(filename)
  P = RtsneWithP::Rtsne(X, perplexity = perplexity, theta = 0, max_iter = 0)$P
  sscore = neMDBD::singularity_score_compute(Y,P)
  sscore_list[[i]] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(1)] = Y[,1]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(2)] = Y[,2]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),3] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),4] = perplexity
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),6] = (sscore > 40000)
}
options(scipen = -1)
my_labeller = function(labels) {return(paste("Perplexity", labels))}

# Create the plot of singularity score (dichotomized)
colors = c('#FF1F5B','#A0B1BA')
plot_mat$bi_label = factor(plot_mat$bi_label, levels = c('1', '0'))
p_dicho = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = as.factor(bi_label), size = as.factor(bi_label))) +
  geom_point(show.legend = TRUE)+
  scale_size_manual(values = c("1" = 1, "0" = 1))+  
  guides(size = "none")  +
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  xlab('tSNE1') + ylab('tSNE2') +
  scale_color_manual(values = colors,name = "Dichotomized\nSingularity\nScore", labels = c('score\n> 40000', 'otherwise')) + 
  guides(color = guide_legend(override.aes = list(size = 2)))

# plot by label
p_label = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = label)) +
  geom_point(size = 1.2, show.legend = TRUE) +
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  xlab('tSNE1') + ylab('tSNE2') +
  labs(color = "Cell Type") + 
  guides(color = guide_legend(override.aes = list(size = 2)))

P = egg::ggarrange(p_label, p_dicho, nrow = 1)
ggsave(P, filename = './plots/SingleCell/brain/brain_plot.png', width = 10, height = 8)

#######################################
##### degree of FI discontinuity decreases in perplexity

perplexity_vec = c(seq(5,150,5))
load('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_5.RData')
sscore_list = vector("list", length = length(perplexity_vec))
plot_mat = data.frame(y1 = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      y2 = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      sscore = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      perplexity = rep(0,length(perplexity_vec)*dim(Y)[1]))
pcol_list = vector('list', length = length(perplexity_vec))
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  load(paste('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_',perplexity,'_Y.RData', sep = ''))
  P = RtsneWithP::Rtsne(X, perplexity = perplexity, theta = 0, max_iter = 0)$P
  sscore = neMDBD::singularity_score_compute(Y,P)#
  sscore_list[[i]] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(1,2)] = Y
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),3] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),4] = perplexity
}

q = 0.05
mean_1_sscore_q = numeric(length(perplexity_vec))
for (i in 1:length(perplexity_vec)){
  sscore = sscore_list[[i]]
  mean_1_sscore_q[[i]] = mean(sscore[sscore>quantile(sscore, 1-q)])
}
p_sscore = ggplot(data = data.frame(mean = mean_1_sscore_q, perplexity = perplexity_vec)) +
  geom_point(aes(x = perplexity, y = mean)) + ggtitle('Degree of FI Discontinuity')+ theme_classic()+
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.5), 
    panel.grid.minor = element_line(color = "grey95", size = 0.25) 
  ) + ylab('Mean of top 5% singularity score') + xlab('Perplexity') + 
  geom_point(data = data.frame(x = c(5,95), y = c(mean_1_sscore_q[c(1,19)])), aes(x = x, y = y), size = 2, colour='#FF1F5B')

#######################################
##### neighborhood preservation increases in perplexity

get_cor_score = function(X, Y, dist_X, nn_mat) {
  dist_Y = as.matrix(dist(Y))
  score_vec = sapply(c(1:dim(X)[1]), function(i) {
    nn_ind = nn_mat[i,]
    return(cor(c(dist_X[i,nn_ind]),c(dist_Y[i,nn_ind])))
  })
  return(score_vec)
}

perplexity_vec = c(seq(5,150,5))
load('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_5.RData')
cor_score_vec = c(rep(0, length(perplexity_vec)))
pca_result = prcomp(X, rank = 50)
dist_X = as.matrix(dist(pca_result$x[,1:50]))
knn_result = FNN::get.knn(pca_result$x[,1:50], round(dim(X)[1]/5))
nn_mat = knn_result$nn.index
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  filename = paste('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_',perplexity,'.RData', sep = '')
  load(filename)
  cor_score_vec[i] = median(get_cor_score(pca_result$x[,1:50], Y, dist_X, nn_mat))
}

p_nnpre = ggplot(data = data.frame(score = cor_score_vec, perplexity = perplexity_vec)) +
  geom_point(aes(x = perplexity, y = score)) + ggtitle('Neighborhood Preservation')+ theme_classic()+
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.5), 
    panel.grid.minor = element_line(color = "grey95", size = 0.25) 
  ) + ylab('Correlation of NN distances') + xlab('Perplexity') + 
  geom_point(data = data.frame(x = c(5,95), y = c(cor_score_vec[c(1,19)])), aes(x = x, y = y), size = 2, colour='#FF1F5B')

P = egg::ggarrange(p_sscore, p_nnpre, ncol = 1)
ggsave(P, filename = './plots/SingleCell/brain/brain_choose_perplexity.png', width = 5, height = 8)

#######################################
##### tests

## Spearson's test
load('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_5.RData')
label = as.factor(ct.X1)
singularity_score = neMDBD::singularity_score_compute(Y, P)
# output mean list for each group
mean_vec = list()
mean_mat = matrix(0, nrow = dim(Y)[1], ncol = dim(Y)[2])
label_unique = c('Astrocytes', 'Endothelial Cells', 'Excitatory Neurons', 'Inhibitory Neurons', 'Microglia', 'Oligodendrocytes')
for (cl in label_unique){
  mean_vec = c(mean_vec, list(colMeans(Y[label == cl,])))
  mean_mat[label == cl,] = t(replicate(sum(label == cl),colMeans(Y[label == cl,])))
}
# output distance to respective group center
dist2center = sqrt(rowSums((Y-mean_mat)**2))
# perform Spearman rank test in R
p_values = numeric()
for (cl in label_unique){
  test_result = cor.test(singularity_score[label == cl], dist2center[label == cl], method = "spearman")
  p_values = c(p_values, test_result$p.value)
}
print('Perplexity 5: p values are (may consider Holm-Bonferroni correction)')
print(p_values)
# [1] 0.21318423 0.55936717 0.74611437 0.04700707
# [5] 0.37351260 0.63637247
print(p.adjust(p_values, method = "holm"))
# [1] 1.0000000 1.0000000 1.0000000 0.2820424
# [5] 1.0000000 1.0000000

load('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_95.RData')
singularity_score = neMDBD::singularity_score_compute(Y, P)
# output mean list for each group
mean_vec = list()
mean_mat = matrix(0, nrow = dim(Y)[1], ncol = dim(Y)[2])
label_unique = c('Astrocytes', 'Endothelial Cells', 'Excitatory Neurons', 'Inhibitory Neurons', 'Microglia', 'Oligodendrocytes')
for (cl in label_unique){
  mean_vec = c(mean_vec, list(colMeans(Y[label == cl,])))
  mean_mat[label == cl,] = t(replicate(sum(label == cl),colMeans(Y[label == cl,])))
}
# output distance to respective group center
dist2center = sqrt(rowSums((Y-mean_mat)**2))
# perform Spearman rank test in R
p_values = numeric()
for (cl in label_unique){
  test_result = cor.test(singularity_score[label == cl], dist2center[label == cl], method = "spearman")
  p_values = c(p_values, test_result$p.value)
}
print('Perplexity 150: p values are (may consider Holm-Bonferroni correction)')
print(p_values)
# [1] 5.047729e-01 6.510493e-02 1.172935e-05
# [4] 6.311108e-02 1.223032e-07 3.888304e-03
print(p.adjust(p_values, method = "holm"))
# [1] 5.047729e-01 1.893332e-01 5.864677e-05
# [4] 1.893332e-01 7.338190e-07 1.555322e-02

## F test for local regression
load('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_5.RData')
singularity_score = neMDBD::singularity_score_compute(Y, P)
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = singularity_score, one = rep(1, dim(Y)[1]))
fit_full = loess(score ~ y1 * y2, span = 0.8, degree = 1, data = data)
fit_reduced = loess(score ~ -y1-y2, span = 1, degree = 0, data = data)
print(anova(fit_full,fit_reduced))
# Model 1: loess(formula = score ~ y1 * y2, data = data, span = 0.8, degree = 1)
# Model 2: loess(formula = score ~ -y1 - y2, data = data, span = 1, degree = 0)
# 
# Analysis of Variance:   denominator df 3611.78
# 
# ENP        RSS F-value Pr(>F)
# [1,] 4.89 8.5837e+13               
# [2,] 1.65 8.6084e+13  1.8661 0.1013

load('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_95.RData')
singularity_score = neMDBD::singularity_score_compute(Y, P)
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = singularity_score, one = rep(1, dim(Y)[1]))
fit_full = loess(score ~ y1 * y2, span = 0.8, degree = 1, data = data)
fit_reduced = loess(score ~ -y1-y2, span = 1, degree = 0, data = data)
print(anova(fit_full,fit_reduced))
# Model 1: loess(formula = score ~ y1 * y2, data = data, span = 0.8, degree = 1)
# Model 2: loess(formula = score ~ -y1 - y2, data = data, span = 1, degree = 0)
# 
# Analysis of Variance:   denominator df 3611.81
# 
# ENP        RSS F-value    Pr(>F)    
# [1,] 4.87 1.1843e+12                      
# [2,] 1.44 1.2290e+12  23.821 < 2.2e-16 ***
#   ---
#   Signif. codes:  
#   0 ‘***’ 1e-03 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Permutation test for local regression
load('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_5.RData')
singularity_score = neMDBD::singularity_score_compute(Y, P)
set.seed(1)
n = 2000
F_vec = numeric(n)
RSS_reduced = (dim(Y)[1]-1) * var(singularity_score)
for (i in c(1:n)){
  data = data.frame(y1 = Y[,1], y2 = Y[,2], score = sample(singularity_score))
  fit_full = loess(score ~ y1 * y2, span = 0.8, degree = 1, data = data)
  RSS_full = sum(fit_full$residuals**2)
  F_vec[i] = (RSS_reduced - RSS_full) / RSS_full
}
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = (singularity_score))
fit_full = loess(score ~ y1 * y2, span = 0.8, degree = 1, data = data)
RSS_full = sum(fit_full$residuals**2)
F_stat = (RSS_reduced - RSS_full) / RSS_full
print(sum(F_stat<F_vec)/length(F_vec))
# [1] 0.0610

load('./data/SingleCell/small_atac_gene_activity_10geno_tsne_perplexity_95.RData')
singularity_score = neMDBD::singularity_score_compute(Y, P)
set.seed(1)
n = 2000
F_vec = numeric(n)
RSS_reduced = (dim(Y)[1]-1) * var(singularity_score)
for (i in c(1:n)){
  data = data.frame(y1 = Y[,1], y2 = Y[,2], score = sample(singularity_score))
  fit_full = loess(score ~ y1 * y2, span = 0.8, degree = 1, data = data)
  RSS_full = sum(fit_full$residuals**2)
  F_vec[i] = (RSS_reduced - RSS_full) / RSS_full
}
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = (singularity_score))
fit_full = loess(score ~ y1 * y2, span = 0.8, degree = 1, data = data)
RSS_full = sum(fit_full$residuals**2)
F_stat = (RSS_reduced - RSS_full) / RSS_full
print(sum(F_stat<F_vec)/length(F_vec))
# [1] 0.0000

########################################
# Mouse mammary epithelial cell
## Updated on May 9, 2025
########################################

#######################################
##### compare small perplexity and large perplexity

perplexity_vec = c(30, 175)
load('./data/SingleCell/BachMammary_pca_50.RData') # The mouse mammary epithelial data is from 'BachMammaryData' in the R package 'scRNAseq'.
load('./data/SingleCell/BachMammary_label.RData')
X = X_pca_50
load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_',perplexity_vec[1],'.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_',perplexity_vec[1],'.RData',sep=''))
sscore_list = vector("list", length = length(perplexity_vec))
plot_mat = data.frame(y1 = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      y2 = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      sscore = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      perplexity = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      label = rep(as.factor(factor(bachmammary_label,labels = c("G1", "G2", "L1", "L2", "NP1", "NP2", "PI1", "PI2"))),length(perplexity_vec)),
                      bi_label = rep(0,length(perplexity_vec)*dim(Y)[1]))
pcol_list = vector('list', length = length(perplexity_vec))
# the following for-loop may take a few minutes to execute
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_',perplexity,'.RData',sep=''))
  load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_',perplexity,'.RData',sep=''))
  # P = RtsneWithP::Rtsne(X, perplexity = perplexity, theta = 0, max_iter = 0)$P
  # sscore = neMDBD::singularity_score_compute(Y,P)
  sscore_list[[i]] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(1)] = Y[,1]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(2)] = Y[,2]
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),3] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),4] = perplexity
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),6] = (sscore > 450000)
}
options(scipen = -1)
my_labeller = function(labels) {return(paste("Perplexity", labels))}

# Create the plot of singularity score (dichotomized)
colors = c('#FF1F5B','#A0B1BA')
plot_mat$bi_label = factor(plot_mat$bi_label, levels = c('1', '0'))
p_dicho = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = as.factor(bi_label), size = as.factor(bi_label))) +
  geom_point(show.legend = TRUE)+
  scale_size_manual(values = c("1" = 1, "0" = 1))+  
  guides(size = "none")  +
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  xlab('tSNE1') + ylab('tSNE2') +
  scale_color_manual(values = colors,name = "Dichotomized\nSingularity\nScore", labels = c('score\n> 450000', 'otherwise')) + 
  guides(color = guide_legend(override.aes = list(size = 2)))

# plot by label
p_label = ggplot(data = data.frame(plot_mat), aes(x = y1, y = y2, color = label)) +
  geom_point(size = 1.2, show.legend = TRUE) +
  facet_wrap(~perplexity, nrow = 2, scales = "free", labeller = as_labeller(my_labeller)) +
  xlab('tSNE1') + ylab('tSNE2') +
  labs(color = "Cell Type") + 
  guides(color = guide_legend(override.aes = list(size = 2)))

P = egg::ggarrange(p_label, p_dicho, nrow = 1)
ggsave(P, filename = './plots/SingleCell/mammary/mammary_plot.png', width = 10, height = 8)

#######################################
##### degree of FI discontinuity decreases in perplexity

perplexity_vec = c(25,30,seq(50,500,25))
load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_5.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_5.RData',sep=''))
sscore_list = vector("list", length = length(perplexity_vec))
plot_mat = data.frame(y1 = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      y2 = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      sscore = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      perplexity = rep(0,length(perplexity_vec)*dim(Y)[1]))
pcol_list = vector('list', length = length(perplexity_vec))
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_',perplexity,'.RData',sep=''))
  load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_',perplexity,'.RData',sep=''))
  # P = RtsneWithP::Rtsne(X, perplexity = perplexity, theta = 0, max_iter = 0)$P
  # sscore = neMDBD::singularity_score_compute(Y,P)#
  sscore_list[[i]] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),c(1,2)] = Y
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),3] = sscore
  plot_mat[((i-1)*dim(X)[1]+1):(i*dim(X)[1]),4] = perplexity
}

q = 0.05
mean_1_sscore_q = numeric(length(perplexity_vec))
for (i in 1:length(perplexity_vec)){
  sscore = sscore_list[[i]]
  mean_1_sscore_q[[i]] = mean(sscore[sscore>quantile(sscore, 1-q)])
}
p_sscore = ggplot(data = data.frame(mean = mean_1_sscore_q, perplexity = perplexity_vec)) +
  geom_point(aes(x = perplexity, y = mean)) + ggtitle('Degree of FI Discontinuity')+ theme_classic()+
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.5), 
    panel.grid.minor = element_line(color = "grey95", size = 0.25) 
  ) + ylab('Mean of top 5% singularity score') + xlab('Perplexity') + 
  geom_point(data = data.frame(x = c(30,175), y = c(mean_1_sscore_q[c(2,8)])), aes(x = x, y = y), size = 2, colour='#FF1F5B')

#######################################
##### neighborhood preservation increases in perplexity

get_cor_score = function(X, Y, dist_X, nn_mat) {
  dist_Y = as.matrix(dist(Y))
  score_vec = sapply(c(1:dim(X)[1]), function(i) {
    nn_ind = nn_mat[i,]
    return(cor(c(dist_X[i,nn_ind]),c(dist_Y[i,nn_ind])))
  })
  return(score_vec)
}

perplexity_vec = c(25,30,seq(50,500,25))
load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_5.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_5.RData',sep=''))
cor_score_vec = c(rep(0, length(perplexity_vec)))
pca_result = X_pca_50
dist_X = as.matrix(dist(pca_result))
knn_result = FNN::get.knn(pca_result, round(dim(X)[1]/5))
nn_mat = knn_result$nn.index
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_',perplexity,'.RData',sep=''))
  load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_',perplexity,'.RData',sep=''))
  cor_score_vec[i] = median(get_cor_score(X_pca_50, Y, dist_X, nn_mat))
  print(i)
}

p_nnpre = ggplot(data = data.frame(score = cor_score_vec, perplexity = perplexity_vec)) +
  geom_point(aes(x = perplexity, y = score)) + ggtitle('Neighborhood Preservation')+ theme_classic()+
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.5), 
    panel.grid.minor = element_line(color = "grey95", size = 0.25) 
  ) + ylab('Correlation of NN distances') + xlab('Perplexity') + 
  geom_point(data = data.frame(x = c(30,175), y = c(cor_score_vec[c(2,8)])), aes(x = x, y = y), size = 2, colour='#FF1F5B')

P = egg::ggarrange(p_sscore, p_nnpre, ncol = 1)
ggsave(P, filename = './plots/SingleCell/mammary/mammary_choose_perplexity.png', width = 5, height = 8)

#######################################
##### tests

## Spearson's test
load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_5.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_5.RData',sep=''))
label = as.factor(factor(bachmammary_label,labels = c("G1", "G2", "L1", "L2", "NP1", "NP2", "PI1", "PI2")))
singularity_score = sscore
# singularity_score = neMDBD::singularity_score_compute(Y, P)
# output mean list for each group
mean_vec = list()
mean_mat = matrix(0, nrow = dim(Y)[1], ncol = dim(Y)[2])
label_unique = unique(label)
for (cl in label_unique){
  mean_vec = c(mean_vec, list(colMeans(Y[label == cl,])))
  mean_mat[label == cl,] = t(replicate(sum(label == cl),colMeans(Y[label == cl,])))
}
# output distance to respective group center
dist2center = sqrt(rowSums((Y-mean_mat)**2))
# perform Spearman rank test in R
p_values = numeric()
for (cl in label_unique){
  test_result = cor.test(singularity_score[label == cl], dist2center[label == cl], method = "spearman")
  p_values = c(p_values, test_result$p.value)
}
print('Perplexity 5: p values are (may consider Holm-Bonferroni correction)')
print(p_values)
# [1] 9.962285e-01 5.565370e-02 6.469455e-01 1.467752e-02 2.739562e-25 1.526189e-08 5.188365e-01 2.184748e-03
print(p.adjust(p_values, method = "holm"))
# [1] 1.000000e+00 2.226148e-01 1.000000e+00 7.338759e-02 2.191649e-24 1.068332e-07 1.000000e+00 1.310849e-02

load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_30.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_30.RData',sep=''))
label = as.factor(factor(bachmammary_label,labels = c("G1", "G2", "L1", "L2", "NP1", "NP2", "PI1", "PI2")))
singularity_score = sscore
# singularity_score = neMDBD::singularity_score_compute(Y, P)
# output mean list for each group
mean_vec = list()
mean_mat = matrix(0, nrow = dim(Y)[1], ncol = dim(Y)[2])
label_unique = unique(label)
for (cl in label_unique){
  mean_vec = c(mean_vec, list(colMeans(Y[label == cl,])))
  mean_mat[label == cl,] = t(replicate(sum(label == cl),colMeans(Y[label == cl,])))
}
# output distance to respective group center
dist2center = sqrt(rowSums((Y-mean_mat)**2))
# perform Spearman rank test in R
p_values = numeric()
for (cl in label_unique){
  test_result = cor.test(singularity_score[label == cl], dist2center[label == cl], method = "spearman")
  p_values = c(p_values, test_result$p.value)
}
print('Perplexity 30: p values are (may consider Holm-Bonferroni correction)')
print(p_values)
# [1] 5.881821e-02 2.149524e-01 1.161762e-02 5.002726e-06 4.400030e-38 1.974548e-13 7.725237e-02 1.809142e-04
print(p.adjust(p_values, method = "holm"))
# [1] 1.764546e-01 2.149524e-01 4.647049e-02 3.001635e-05 3.520024e-37 1.382184e-12 1.764546e-01 9.045708e-04

load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_175.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_175.RData',sep=''))
label = as.factor(factor(bachmammary_label,labels = c("G1", "G2", "L1", "L2", "NP1", "NP2", "PI1", "PI2")))
singularity_score = sscore
# singularity_score = neMDBD::singularity_score_compute(Y, P)
# output mean list for each group
mean_vec = list()
mean_mat = matrix(0, nrow = dim(Y)[1], ncol = dim(Y)[2])
label_unique = unique(label)
for (cl in label_unique){
  mean_vec = c(mean_vec, list(colMeans(Y[label == cl,])))
  mean_mat[label == cl,] = t(replicate(sum(label == cl),colMeans(Y[label == cl,])))
}
# output distance to respective group center
dist2center = sqrt(rowSums((Y-mean_mat)**2))
# perform Spearman rank test in R
p_values = numeric()
for (cl in label_unique){
  test_result = cor.test(singularity_score[label == cl], dist2center[label == cl], method = "spearman")
  p_values = c(p_values, test_result$p.value)
}
print('Perplexity 175: p values are (may consider Holm-Bonferroni correction)')
print(p_values)
# [1]  3.623494e-01  4.702473e-01  7.067602e-14  1.135859e-15 3.278601e-150  5.210038e-76  8.745915e-01  4.573384e-11
print(p.adjust(p_values, method = "holm"))
# [1]  1.000000e+00  1.000000e+00  3.533801e-13  6.815151e-15 2.622881e-149  3.647026e-75  1.000000e+00  1.829354e-10

## F test for local regression
load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_5.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_5.RData',sep=''))
# singularity_score = neMDBD::singularity_score_compute(Y, P)
singularity_score = sscore
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = singularity_score, one = rep(1, dim(Y)[1]))
fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
fit_reduced = loess(score ~ -y1-y2, span = 1, degree = 0, data = data)
print(anova(fit_full,fit_reduced))
# Model 1: loess(formula = score ~ y1 * y2, data = data, span = 0.5, degree = 1)
# Model 2: loess(formula = score ~ -y1 - y2, data = data, span = 1, degree = 0)
# 
# Analysis of Variance:   denominator df 25796.39
# 
# ENP        RSS F-value Pr(>F)
# [1,] 7.10 3.2003e+16               
# [2,] 1.72 3.2010e+16  0.6414 0.7451

load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_30.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_30.RData',sep=''))
# singularity_score = neMDBD::singularity_score_compute(Y, P)
singularity_score = sscore
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = singularity_score, one = rep(1, dim(Y)[1]))
fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
fit_reduced = loess(score ~ -y1-y2, span = 1, degree = 0, data = data)
print(anova(fit_full,fit_reduced))
# Model 1: loess(formula = score ~ y1 * y2, data = data, span = 0.5, degree = 1)
# Model 2: loess(formula = score ~ -y1 - y2, data = data, span = 1, degree = 0)
# 
# Analysis of Variance:   denominator df 25796.92
# 
# ENP        RSS F-value Pr(>F)
# [1,] 6.75 9.0139e+17               
# [2,] 1.79 9.0175e+17  1.1978 0.2979

load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_175.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_175.RData',sep=''))
# singularity_score = neMDBD::singularity_score_compute(Y, P)
singularity_score = sscore
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = singularity_score, one = rep(1, dim(Y)[1]))
fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
fit_reduced = loess(score ~ -y1-y2, span = 1, degree = 0, data = data)
print(anova(fit_full,fit_reduced))
# Model 1: loess(formula = score ~ y1 * y2, data = data, span = 0.5, degree = 1)
# Model 2: loess(formula = score ~ -y1 - y2, data = data, span = 1, degree = 0)
# 
# Analysis of Variance:   denominator df 25797.82
# 
# ENP        RSS F-value    Pr(>F)    
# [1,] 6.17 2.7528e+15                      
# [2,] 2.03 2.7621e+15  11.902 5.381e-14 ***
#   ---
#   Signif. codes: 0 ‘***’ 1e-03 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Permutation test for local regression
load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_5.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_5.RData',sep=''))
# singularity_score = neMDBD::singularity_score_compute(Y, P)
singularity_score = sscore
n = 10000
RSS_reduced = (dim(Y)[1]-1) * var(singularity_score)
# The following code was run on server due to computation time.
# F_vec = numeric(n)
# for (i in c(1:n)){
#   data = data.frame(y1 = Y[,1], y2 = Y[,2], score = sample(singularity_score))
#   fit_full = loess(score ~ y1 * y2, span = 0.8, degree = 1, data = data)
#   RSS_full = sum(fit_full$residuals**2)
#   F_vec[i] = (RSS_reduced - RSS_full) / RSS_full
# }
F_vec <- c(read.table("./data/SingleCell/BachMammary_F_vec_5.txt", header = FALSE, sep = "\t")[,1])
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = (singularity_score))
fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
RSS_full = sum(fit_full$residuals**2)
F_stat = (RSS_reduced - RSS_full) / RSS_full
print(sum(F_stat<F_vec)/length(F_vec))
# [1] 0.8516

load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_30.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_30.RData',sep=''))
# singularity_score = neMDBD::singularity_score_compute(Y, P)
singularity_score = sscore
n = 10000
RSS_reduced = (dim(Y)[1]-1) * var(singularity_score)
# The following code was run on server due to computation time.
# F_vec = numeric(n)
# for (i in c(1:n)){
#   data = data.frame(y1 = Y[,1], y2 = Y[,2], score = sample(singularity_score))
#   fit_full = loess(score ~ y1 * y2, span = 0.8, degree = 1, data = data)
#   RSS_full = sum(fit_full$residuals**2)
#   F_vec[i] = (RSS_reduced - RSS_full) / RSS_full
# }
F_vec <- c(read.table("./data/SingleCell/BachMammary_F_vec_30.txt", header = FALSE, sep = "\t")[,1])
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = (singularity_score))
fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
RSS_full = sum(fit_full$residuals**2)
F_stat = (RSS_reduced - RSS_full) / RSS_full
print(sum(F_stat<F_vec)/length(F_vec))
# [1] 0.2379

load(paste('./data/SingleCell/BachMammaryData_sscore_tSNE_175.RData',sep=''))
load(paste('./data/SingleCell/BachMammaryData_tSNE_perplexity_175.RData',sep=''))
# singularity_score = neMDBD::singularity_score_compute(Y, P)
singularity_score = sscore
n = 10000
RSS_reduced = (dim(Y)[1]-1) * var(singularity_score)
# The following code was run on server due to computation time.
# F_vec = numeric(n)
# for (i in c(1:n)){
#   data = data.frame(y1 = Y[,1], y2 = Y[,2], score = sample(singularity_score))
#   fit_full = loess(score ~ y1 * y2, span = 0.8, degree = 1, data = data)
#   RSS_full = sum(fit_full$residuals**2)
#   F_vec[i] = (RSS_reduced - RSS_full) / RSS_full
# }
F_vec <- c(read.table("./data/SingleCell/BachMammary_F_vec_175.txt", header = FALSE, sep = "\t")[,1])
data = data.frame(y1 = Y[,1], y2 = Y[,2], score = (singularity_score))
fit_full = loess(score ~ y1 * y2, span = 0.5, degree = 1, data = data)
RSS_full = sum(fit_full$residuals**2)
F_stat = (RSS_reduced - RSS_full) / RSS_full
print(sum(F_stat<F_vec)/length(F_vec))
# [1] 0.0000

