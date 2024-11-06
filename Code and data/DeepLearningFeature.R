########################################
# Oct 31, 2024, Zhexuan Liu #############
# Deep learning feature data ###########
########################################

########################################
# CIFAR-10
########################################

library(neMDBD)
library(RtsneWithP)
library(ggplot2)
library(latex2exp)
library(viridis)
library(ROCR)

# We subsampled the test dataset to 5000 points and stored the embedding points (indices stored in 'cifar10_5000_ind.RData').
load('./data/DeepLearningFeatures/embedding/cifar10_perplexity_125.RData')
load('./data/DeepLearningFeatures/embedding/cifar10_5000_label.RData')
load('./data/DeepLearningFeatures/embedding/cifar10_5000_ind.RData')

# # calculate the perturbation score (done in HPC clusters with parallel computing)
# perturbation_score = perturbation_score_compute(X, tsne_out$Y, perplexity, length = 2, approx = 0)
# save(perturbation_score, file = './data/DeepLearningFeatures/perturbation_scores/CIFAR10/len2/perturbation_score.RData')
load('./data/DeepLearningFeatures/perturbation_scores/CIFAR10/len2/perturbation_score.RData')

# plot the labels
p_label = ggplot(data = data.frame(tSNE1 = Y[,1], tSNE2 = Y[,2], label = label)) + 
  geom_point(aes(x = tSNE1, y = tSNE2, color = label), size = 0.2)

# plot the perturbation scores
custom_labels = function(x) {ifelse(x >=30, TeX('$\\geq 30$'), format(x))}
perturbation_score_adjusted = perturbation_score
perturbation_score_adjusted[perturbation_score >= 30] = 30
p_pscore = ggplot(data.frame(tSNE1 = Y[,1], tSNE2 = Y[,2], color = perturbation_score_adjusted), aes(x = tSNE1, y = tSNE2, color = color)) +
  geom_point(size = 0.2) +
  scale_color_viridis_c(name = "Perturbation\nScore", labels = custom_labels)

# plot the entropies of predicted probabilities
library(reticulate)
np = import("numpy")
prob_all = np$load('./data/DeepLearningFeatures/features/CIFAR10/CIFAR10_ResNet18_ce_pretrain_probs_test.npy')
probs = prob_all[ind,]
entropies = sapply(1:dim(probs)[1], function(i) {-sum(probs[i,] * log(probs[i,]))})
custom_labels = function(x) {ifelse(x >=1, TeX('$\\geq$1.0'), format(x))}
entropies_adjusted = entropies
entropies_adjusted[entropies>=1] = 1
p_entropy = ggplot(data.frame(tSNE1 = Y[,1], tSNE2 = Y[,2], color = entropies_adjusted), aes(x = tSNE1, y = tSNE2, color = color)) +
  geom_point(size = 0.2) +
  scale_color_viridis_c(name = "Entropy", labels = custom_labels)

ggsave(p_label, file = './plots/DeepLearningFeatures/cifar10_5000_embedding_label.png', width = 5, height = 4)
ggsave(p_pscore, file = './plots/DeepLearningFeatures/cifar10_5000_embedding_pscore.png', width = 5, height = 4)
ggsave(p_entropy, file = './plots/DeepLearningFeatures/cifar10_5000_embedding_entropy.png', width = 5, height = 4)


########################################
# CIFAR-10 & DTD (OOD detection)
########################################

np = import("numpy")
X_test = np$load('./data/DeepLearningFeatures/features/CIFAR10/CIFAR10_ResNet18_ce_pretrain_features_test.npy')
X_test = as.matrix(X_test)
label_test = np$load('./data/DeepLearningFeatures/features/CIFAR10/CIFAR10_ResNet18_ce_pretrain_labels_test.npy')
label_test = as.factor(label_test)
label_text = c('airplane', 'automobile', 'bird', 'cat', 'deer', 'dog', 'frog', 'horse', 'ship', 'truck')
label = as.numeric(label_test)
for (i in 1:10) {label[label == (i)] = label_text[i]}
label_test = as.factor(label)
X_dtd = np$load('./data/DeepLearningFeatures/features/DTD/CIFAR10_ResNet18_ce_pretrain_features_DTD.npy')
X_dtd = as.matrix(X_dtd)
# subsample each dataset for computation
set.seed(2)
major_ind = sample(c(1:dim(X_test)[1]), 2000, replace = FALSE)
ood_ind = sample(c(1:dim(X_dtd)[1]), 1000, replace = FALSE)
# combine dataset
X = rbind(X_test[major_ind,], X_dtd[ood_ind,])
label = as.factor(c(as.character(label_test[major_ind]), c(rep('DTD (OOD)', length(ood_ind)))))
label = factor(label, levels = c('airplane', 'automobile', 'bird', 'cat', 'deer', 'dog', 'frog', 'horse', 'ship', 'truck', 'DTD (OOD)'))
label_binary = as.factor(c(rep('CIFAR10',length(label_test[major_ind])), c(rep('OOD', length(ood_ind)))))

# plot t-SNE embedding points
load("./data/DeepLearningFeatures/embedding/cifar10_DTD_perplexity_100.RData")
lighter_colors = scales::brewer_pal(palette = "Set3")(12)
lighter_colors = lighter_colors[-c(4,6)]
set.seed(4)
colors = c(sample(lighter_colors), '#FF1F5B')
p_label = ggplot() +  
  geom_point(data = data.frame(tSNE1 = Y[,1], tSNE2 = Y[,2], label = label), aes(x=tSNE1, y=tSNE2, color=label), size = 0.2) +
  scale_color_manual(values = colors,name = "Label")

# calculate perturbation score
# perturbation_score = perturbation_score_compute(X, Y, perplexity = 100, length = 2, approx = 0)
load('./data/DeepLearningFeatures/perturbation_scores/CIFAR10andDTD/len2/perturbation_score.RData')

# plot t-SNE embedding points, perturbation scores, ROC (zoom in)
plot_score = function(X, score, point_size = 0.8, direction = 1){
  plot_df = data.frame(x = X[,1], y = X[,2], score = score)
  p = ggplot() + 
    geom_point(data = data.frame(tSNE1 = X[,1], tSNE2 = X[,2], score = score), aes(x = tSNE1, y = tSNE2, color = score), size = point_size) + 
    scale_color_viridis(direction = direction)
  return(p)
}
plot_zoomin_label_score = function(Y, label, score, id) {
  perturb_score_plot = score
  perturb_score_plot[!id] = mean(perturb_score_plot[id])
  perturb_score_plot[perturb_score_plot > 20] = 20
  custom_labels = function(x) {ifelse(x >= 20, TeX('$\\geq$20'), format(x))}
  P_p = plot_score(Y, perturb_score_plot, point_size = 1.5) +
    labs(x = NULL, y = NULL)
  lighter_colors = scales::brewer_pal(palette = "Set3")(12)
  lighter_colors = lighter_colors[-c(4, 6)]
  set.seed(4)
  colors = c(sample(lighter_colors), '#FF1F5B')
  P_c = ggplot() +
    geom_point(data = data.frame(x = Y[, 1], y = Y[, 2], label = label),
               aes(x = x, y = y, color = label), size = 1.5) +
    scale_color_manual(values = colors, name = "Label") +
    labs(x = NULL, y = NULL) + theme(legend.position = 'None')
  
  p_p = P_p +
    scale_color_viridis(name = "Perturbation\nScore", labels = custom_labels, option = "D", begin = 0, end = 0.95) +
    coord_cartesian(xlim = c(min(Y[id, 1]) - .5, max(Y[id, 1]) + .5),
                    ylim = c(min(Y[id, 2]) - .5, max(Y[id, 2]) + .5))
  p_c = P_c + coord_cartesian(xlim = c(min(Y[id, 1]) - .5, max(Y[id, 1]) + .5),
                              ylim = c(min(Y[id, 2]) - .5, max(Y[id, 2]) + .5)) +
    theme(legend.position = 'None')
  
  perturb_score_zoom = perturb_score[id]
  OOD_id = (label[id] == 'DTD (OOD)')
  pred = prediction(perturb_score_zoom, OOD_id)
  perf = performance(pred, "tpr", "fpr")
  perf.auc = performance(pred, measure = "auc")
  p_roc = ggplot(data = data.frame(x = perf@x.values[[1]], y = perf@y.values[[1]]),aes(x = x, y = y)) + geom_line() +
    theme(plot.title = element_text(size = 10)) +
    geom_text(x = 0.75, y = 0.25, label = paste('AUC =', round(perf.auc@y.values[[1]],4)), hjust = 0.75, vjust = -1, size = 5)+ labs(x = NULL, y = NULL)
  
  P = egg::ggarrange(p_c, p_p, p_roc, nrow = 1)
  return(P)
}

id1 = (Y[,1]>(25) & Y[,1]<(35) & Y[,2]>(-40) & Y[,2]<(-20))
p1 = plot_zoomin_label_score(Y, label, perturbation_score, id1)
id2 = (Y[,1]>(-25) & Y[,1]<(0) & Y[,2]>(22) & Y[,2]<(40))
p2 = plot_zoomin_label_score(Y, label, perturbation_score, id2)
id3 = (Y[,1]>(37) & Y[,1]<(50) & Y[,2]>(-30) & Y[,2]<(-10))
p3 = plot_zoomin_label_score(Y, label, perturbation_score, id3)

ggsave(p_label, file = './plots/DeepLearningFeatures/cifar10_DTD_embedding.png', width = 5, height = 4)
ggsave(p1, file = './plots/DeepLearningFeatures/cifar10_DTD_zoom1.png', width = 15, height = 4)
ggsave(p2, file = './plots/DeepLearningFeatures/cifar10_DTD_zoom2.png', width = 15, height = 4)
ggsave(p3, file = './plots/DeepLearningFeatures/cifar10_DTD_zoom3.png', width = 15, height = 4)
