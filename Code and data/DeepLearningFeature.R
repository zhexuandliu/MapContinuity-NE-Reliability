########################################
# Oct 31, 2024, Zhexuan Liu #############
# Deep learning feature data ###########
########################################

########################################
# CIFAR-10
########################################

########################################
# Perturbation score

library(neMDBD)
library(RtsneWithP)
library(ggplot2)
library(latex2exp)
library(viridis)
library(ROCR)

# We subsampled the test dataset to 5000 points and stored the embedding points (indices stored in 'cifar10_5000_ind.RData').
# set.seed(1)
# ind = sample(1:10000, 5000, replace = FALSE)
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

########################################
# Singularity score

load(paste('./data/DeepLearningFeatures/singularity_scores/cifar10_perplexity_',100,'.RData', sep = ''))
eigen_score_list = vector("list", length = length(perplexity_vec))
plot_mat = data.frame(y1 = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      y2 = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      score_eigen = rep(0,length(perplexity_vec)*dim(Y)[1]),
                      perplexity = rep(0,length(perplexity_vec)*dim(Y)[1]))
pcol_list = vector('list', length = length(perplexity_vec))
for (i in 1:length(perplexity_vec)){
  perplexity = perplexity_vec[i]
  filename = paste('./data/DeepLearningFeatures/singularity_scores/cifar10_perplexity_',perplexity,'.RData', sep = '')
  load(filename)
  # singularity_score = neMDBD::singularity_score_compute(Y,P)
  eigen_score_list[[i]] = singularity_score
  plot_mat[((i-1)*dim(Y)[1]+1):(i*dim(Y)[1]),c(1,2)] = Y
  plot_mat[((i-1)*dim(Y)[1]+1):(i*dim(Y)[1]),3] = singularity_score
  plot_mat[((i-1)*dim(Y)[1]+1):(i*dim(Y)[1]),4] = perplexity
}

q = 0.05
mean_1_eigen_all = numeric(length(perplexity_vec))
mean_1_eigen_q = numeric(length(perplexity_vec))
for (i in 1:length(perplexity_vec)){
  score_eigen = eigen_score_list[[i]]
  mean_1_eigen_all[[i]] = mean(score_eigen)
  mean_1_eigen_q[[i]] = mean(score_eigen[score_eigen>quantile(score_eigen, 1-q)])}
p_dot = ggplot(data = data.frame(y = mean_1_eigen_q, x = perplexity_vec), aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(se = FALSE, size = 0.75, span = 0.5)+ # Plot the points
  ggtitle('Degree of FI Discontinuity')+ theme_classic()+
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.5), 
    panel.grid.minor = element_line(color = "grey95", size = 0.25) 
  ) + ylab('Mean of top 5% singularity score') + xlab('Perplexity') + geom_point(data = data.frame(x = c(5,25,100), y = c(mean_1_eigen_q[c(1,2,5)])),aes(x = x, y = y), pch=0,size=5,colour="blue")  + geom_vline(xintercept = 100,linetype="dotted", color = "red", size=0.8)+ 
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
    strip.text = element_text(size = 16)
  )
ggsave(p_dot, file = './plots/DeepLearningFeatures/degree_of_FI_dicontinuity.png', width = 13.5/7*10/4, height = 3/7*10)


my_labeller <- function(labels) {
  return(paste("Perplexity", labels))
}

plot_mat$hh = 1
id = (plot_mat[,4] ==5 | plot_mat[,4] ==5|plot_mat[,4] ==25|plot_mat[,4] ==100)
id2 = (plot_mat[,3] > 50000) & id
plot_mat$score_eigen_adj = plot_mat$score_eigen
plot_mat$score_eigen_adj[plot_mat$score_eigen>300000] = 300000
custom_labels <- function(x) {
  ifelse(x >=300000, TeX('$\\geq$3e+05'), format(x))
}
p <- ggplot(data = plot_mat[id,], aes(x = y1, y = y2)) +
  geom_point(data = subset(plot_mat, id2), aes(x = y1, y = y2, shape = "Special"), size = 1.25, colour="red") +
  geom_point(aes(color = (score_eigen_adj)), size = 0.4) +  # Regular points  # Red circles
  scale_shape_manual(values = c(Special = 1), labels = 'with severe\n FI discontinuity') +  # Define shape for "special" points, no legend title
  scale_color_viridis(direction = 1, trans = 'log10', name = "Singularity\nScore",
                      labels = custom_labels) +  # Color scale for eigenscore
  facet_wrap(perplexity ~ ., nrow = 1, scales = "free", labeller = as_labeller(my_labeller)) + 
  xlab('tSNE1') + ylab('tSNE2') +
  guides(shape = guide_legend(title =TeX("Score$>5\\times 10^4$")))+ 
  theme_minimal()+
  theme(strip.text = element_text(colour = "blue", face = "bold.italic", size = 16),
        strip.background = element_blank())+ 
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16)
  )
ggsave(p, file = './plots/DeepLearningFeatures/compare_different_perplexity.png', width = 14.3/7*10/4*3, height = 3/7*10)

########################################
# Compare with Kernel PCA and One-class SVM
# Updated on May 10, 2025

test_OOD_KPCA_score = np$load("./data/DeepLearningFeatures/KPCA/KPCA_CoP_In.npy")
dtd_OOD_KPCA_score = np$load("./data/DeepLearningFeatures/KPCA/KPCA_CoP_OOD.npy")
KPCA_score = c(test_OOD_KPCA_score, dtd_OOD_KPCA_score)
test_OOD_OCSVM_score = np$load("./data/DeepLearningFeatures/OCSVM/OCSVM_In.npy")
dtd_OOD_OCSVM_score = np$load("./data/DeepLearningFeatures/OCSVM/OCSVM_OOD.npy")
OCSVM_score = c(test_OOD_OCSVM_score, dtd_OOD_OCSVM_score)
label = c(rep("InD", 2000), rep('DTD (OOD)', 1000))

KPCA_score_zoom = KPCA_score[id1]
OOD_id = (label[id1] == 'DTD (OOD)')
pred = prediction(KPCA_score_zoom, OOD_id)
perf = performance(pred, "tpr", "fpr")
perf.auc = performance(pred, measure = "auc")
ROC_curve_KPCA = ggplot(data = data.frame(x = perf@x.values[[1]], y = perf@y.values[[1]]),aes(x = x, y = y)) + geom_line() +
  theme(plot.title = element_text(size = 10)) +
  geom_text(x = 0.75, y = 0.25,
            label = paste('AUC =', round(perf.auc@y.values[[1]],4)),
            hjust = 0.75, vjust = -1, size = 5)+
  labs(x = NULL, y = NULL)+
  theme_minimal()+
  theme_minimal() +
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16)
  )  + xlab('False Positive Rate') + ylab('True Positive Rate')
OCSVM_score_zoom = -OCSVM_score[id1]
OOD_id = (label[id1] == 'DTD (OOD)')
pred = prediction(OCSVM_score_zoom, OOD_id)
perf = performance(pred, "tpr", "fpr")
perf.auc = performance(pred, measure = "auc")
ROC_curve_OCSVM = ggplot(data = data.frame(x = perf@x.values[[1]], y = perf@y.values[[1]]),aes(x = x, y = y)) + geom_line() +
  theme(plot.title = element_text(size = 10)) +
  geom_text(x = 0.75, y = 0.25,
            label = paste('AUC =', round(perf.auc@y.values[[1]],4)),
            hjust = 0.75, vjust = -1, size = 5)+
  labs(x = NULL, y = NULL)+
  theme_minimal()+
  theme_minimal() +
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16)
  )  + xlab('False Positive Rate') + ylab('True Positive Rate')
P = egg::ggarrange(ROC_curve_KPCA + ggtitle("Kernel PCA"), 
                   ROC_curve_OCSVM + ggtitle("One-class SVM"), ncol = 2)
ggsave(P, file = "./plots/DeepLearningFeatures/compare_1.png", width = 13.5/7*10/2, height = 3.6/7*10)

KPCA_score_zoom = KPCA_score[id2]
OOD_id = (label[id2] == 'DTD (OOD)')
pred = prediction(KPCA_score_zoom, OOD_id)
perf = performance(pred, "tpr", "fpr")
perf.auc = performance(pred, measure = "auc")
ROC_curve_KPCA = ggplot(data = data.frame(x = perf@x.values[[1]], y = perf@y.values[[1]]),aes(x = x, y = y)) + geom_line() +
  theme(plot.title = element_text(size = 10)) +
  geom_text(x = 0.75, y = 0.25,
            label = paste('AUC =', round(perf.auc@y.values[[1]],4)),
            hjust = 0.75, vjust = -1, size = 5)+
  labs(x = NULL, y = NULL)+
  theme_minimal()+
  theme_minimal() +
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16)
  )  + xlab('False Positive Rate') + ylab('True Positive Rate')
OCSVM_score_zoom = -OCSVM_score[id2]
OOD_id = (label[id2] == 'DTD (OOD)')
pred = prediction(OCSVM_score_zoom, OOD_id)
perf = performance(pred, "tpr", "fpr")
perf.auc = performance(pred, measure = "auc")
ROC_curve_OCSVM = ggplot(data = data.frame(x = perf@x.values[[1]], y = perf@y.values[[1]]),aes(x = x, y = y)) + geom_line() +
  theme(plot.title = element_text(size = 10)) +
  geom_text(x = 0.75, y = 0.25,
            label = paste('AUC =', round(perf.auc@y.values[[1]],4)),
            hjust = 0.75, vjust = -1, size = 5)+
  labs(x = NULL, y = NULL)+
  theme_minimal()+
  theme_minimal() +
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16)
  )  + xlab('False Positive Rate') + ylab('True Positive Rate')
P = egg::ggarrange(ROC_curve_KPCA + ggtitle("Kernel PCA"), 
                   ROC_curve_OCSVM + ggtitle("One-class SVM"), ncol = 2)
ggsave(P, file = "./plots/DeepLearningFeatures/compare_2.png", width = 13.5/7*10/2, height = 3.6/7*10)

KPCA_score_zoom = KPCA_score[id3]
OOD_id = (label[id3] == 'DTD (OOD)')
pred = prediction(KPCA_score_zoom, OOD_id)
perf = performance(pred, "tpr", "fpr")
perf.auc = performance(pred, measure = "auc")
ROC_curve_KPCA = ggplot(data = data.frame(x = perf@x.values[[1]], y = perf@y.values[[1]]),aes(x = x, y = y)) + geom_line() +
  theme(plot.title = element_text(size = 10)) +
  geom_text(x = 0.75, y = 0.25,
            label = paste('AUC =', round(perf.auc@y.values[[1]],4)),
            hjust = 0.75, vjust = -1, size = 5)+
  labs(x = NULL, y = NULL)+
  theme_minimal()+
  theme_minimal() +
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16)
  )  + xlab('False Positive Rate') + ylab('True Positive Rate')
OCSVM_score_zoom = -OCSVM_score[id3]
OOD_id = (label[id3] == 'DTD (OOD)')
pred = prediction(OCSVM_score_zoom, OOD_id)
perf = performance(pred, "tpr", "fpr")
perf.auc = performance(pred, measure = "auc")
ROC_curve_OCSVM = ggplot(data = data.frame(x = perf@x.values[[1]], y = perf@y.values[[1]]),aes(x = x, y = y)) + geom_line() +
  theme(plot.title = element_text(size = 10)) +
  geom_text(x = 0.75, y = 0.25,
            label = paste('AUC =', round(perf.auc@y.values[[1]],4)),
            hjust = 0.75, vjust = -1, size = 5)+
  labs(x = NULL, y = NULL)+
  theme_minimal()+
  theme_minimal() +
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16)
  )  + xlab('False Positive Rate') + ylab('True Positive Rate')
P = egg::ggarrange(ROC_curve_KPCA + ggtitle("Kernel PCA"), 
                   ROC_curve_OCSVM + ggtitle("One-class SVM"), ncol = 2)
ggsave(P, file = "./plots/DeepLearningFeatures/compare_3.png", width = 13.5/7*10/2, height = 3.6/7*10)