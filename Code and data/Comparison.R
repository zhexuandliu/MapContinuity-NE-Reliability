########################################
# May 9, 2025, Zhexuan Liu #############
# Comparison with other methods ########
########################################

library(ggplot2)
library(viridis)
library(neMDBD)
library(Rtsne)

########################################
# Comparison with perturbation score
########################################

# Data generation
dataset_name = 'swissroll'
dim = 3
generate_swiss_roll = function(n_samples = 1500, noise = 0.05) {
  t = runif(n_samples, 0, 5*pi)
  x = t * cos(t)
  y = t * sin(t)
  z = runif(n_samples, -1, 1) * 1 
  x = x + rnorm(n_samples, mean = 0, sd = noise)
  y = y + rnorm(n_samples, mean = 0, sd = noise)
  z = z + rnorm(n_samples, mean = 0, sd = noise)
  list(x = x, y = y, z = z, t = t)
}
set.seed(5)
X = generate_swiss_roll(n_sample = 1000)
t = X$t
X = cbind(X$x, X$y, X$z) 
PCA_x = prcomp(X) 
colnames(X) = NULL
rownames(X) = NULL
Y = Rtsne(X, perplexity = 150, theta = 0, max_iter = 5000, Y_init = PCA_x$x[,1:2])
Y = Y$Y

# scDEED
library(scDEED)
library(Seurat)
library(gridExtra)
library(dplyr)
library(patchwork)
library(VGAM)
library(gplots)
library(ggplot2)
library(pracma)
library(resample)
library(foreach)
library(distances)
library(utils)
library(doParallel)
data = X
rownames(data) = paste0("Cell", 1:dim(X)[1])
colnames(data) = paste0("Gene", 1:dim(X)[2])
data = CreateSeuratObject(counts = t(X))
data = NormalizeData(data)
data = ScaleData(data)
data = FindVariableFeatures(data, selection.method = 'dispersion')
data = SetAssayData(data, slot = "scale.data", new.data = t(X))
data = RunPCA(data, npcs = 3)
perp=150
data = RunTSNE(data, perplexity = perp)
tsne_embeddings = Embeddings(data, reduction = "tsne")
modified_tsne_embeddings = Y
rownames(modified_tsne_embeddings) = rownames(tsne_embeddings)
colnames(modified_tsne_embeddings) = colnames(tsne_embeddings)
data[["tsne"]]@cell.embeddings = modified_tsne_embeddings
result = scDEED(data, K = 2, reduction.method = 'tsne', perplexity = perp, rerun = F, dubious_cutoff = 0.05, trustworthy_cutoff = 0.95)
dubious_cells = result$full_results$dubious_cells[result$full_results$perplexity==as.character(perp)]
dubious_cells = as.numeric(strsplit(dubious_cells, ',')[[1]])
trustworthy_cells =  result$full_results$trustworthy_cells[result$full_results$perplexity==as.character(perp)]
trustworthy_cells = as.numeric(strsplit(trustworthy_cells, ',')[[1]])
scDEED_output = c(rep('unselected', 1000))
scDEED_output[dubious_cells] = 'dubious'
scDEED_output[trustworthy_cells] = 'trustworthy'
scDEED_output = as.factor(scDEED_output)
p_scDEED = ggplot(data = data.frame(x = Y[,1], y = Y[,2], t = scDEED_output)) + geom_point(aes(x = x, y = y, color = t), size = 1.2)+ scale_color_manual(values = c('#FF1F5B','#009ADE','gray')) + xlab('tSNE1') + ylab('tSNE2')+ 
  theme_minimal()+
  theme(
    text = element_text(size = 16),  
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 16), 
    plot.title = element_text(size = 20),  
    plot.subtitle = element_text(size = 20),  
    plot.caption = element_text(size = 14),  
    strip.text = element_text(size = 16)
  )+ labs(color = "Label")+ ggtitle('t-SNE Embedding with Labeling by scDEED')


# EMBEDR
EMBEDR_Pvalues = c(read.table("./data/Comparison/EMBEDR_PValues.txt", header = FALSE, sep = "\t"))$V1
df = data.frame(x = Y[,1], y = Y[,2], color_value = EMBEDR_Pvalues)
df$color_range = cut(log10(df$color_value), 
                      breaks = c(-5, -3, -2, -1, 0), 
                      labels = c("1e-5 to 1e-3", "1e-3 to 1e-2", "1e-2 to 1e-1", "1e-1 to 1"),
                      include.lowest = TRUE)
range_colors = c("1e-5 to 1e-3" = "green", 
                  "1e-3 to 1e-2" = "orange", 
                  "1e-2 to 1e-1" = "blue", 
                  "1e-1 to 1" = "purple")
p_EMBEDR = ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color = log10(color_value)), size = 1.2) +
  scale_color_gradientn(colors = c(rgb(77/255, 154/255, 117/255), 
                                   rgb(199/255, 229/255, 221/255), 
                                   rgb(238/255, 214/255, 190/255), 
                                   rgb(193/255, 105/255, 29/255), 
                                   rgb(198/255, 223/255, 239/255), 
                                   rgb(61/255, 113/255, 175/255), 
                                   rgb(237/255, 222/255, 236/255), 
                                   rgb(190/255, 127/255, 185/255)),
                        values = scales::rescale(c(-5, -3.00000001, -3, -2.00000001, -2, -1.00000001, -1, 0)),
                        breaks = c(-5, -3, -2, -1, 0),
                        labels = c("1.0e-05", "1.0e-03", "1.0e-02", "1.0e-01", "1.0e+00"),
                        limits = c(-5, 0),
                        oob = scales::squish) +
  theme_minimal() +
  theme(legend.key.height = unit(1.5, "cm")) +
  labs(color = "EMBEDR\np-value") + 
  ggtitle('t-SNE Embedding with p-values by EMBEDR')+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 16), 
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16)
  )+ xlab('tSNE1') + ylab('tSNE2')

# DynamicViz
variance_scores = c(read.table("./data/Comparison/variance_scores.txt", header = FALSE, sep = "\t"))$V1
df = data.frame(x = Y[,1], y = Y[,2], color_value = variance_scores)
p_DynamicViz = ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color = color_value), size = 1.2) +
  scale_color_gradientn(colors = c("blue", "lightblue", "white", "lightpink", "red"),
                        values = scales::rescale(c(0.01, 0.05, 0.08, 0.1, 0.12)/0.12*0.132),
                        limits = c(0, 0.132))+
  theme_minimal() +
  theme(legend.key.height = unit(1.5, "cm")) +
  labs(color = "variance\nscore") + 
  ggtitle('t-SNE Embedding with scores by DynamicViz')+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 16), 
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16)
  )+ xlab('tSNE1') + ylab('tSNE2')

# Perturbation Score
load("./data/Comparison/perturbation_score.RData")
plot_df = data.frame(x = Y[,1], y = Y[,2], score = perturbation_score)
plot_df$point_size = ifelse(perturbation_score > 8, 3, 1.2)
p_pscore = ggplot() + 
  geom_point(data = plot_df, aes(x = x, y = y, color = score, size = point_size), 
             show.legend = c(color = TRUE, size = FALSE)) + xlab('tSNE1') + ylab('tSNE2')+
  scale_size_identity()+
  viridis::scale_color_viridis(
    direction = 1, 
    name = "Perturbation\nScore",
    limits = c(0,18)) +
  theme_minimal()+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 16), 
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16)
  )+ ggtitle('t-SNE Embedding with Perturbation Score')

P = egg::ggarrange(p_pscore,
                   p_EMBEDR+ 
                     theme(plot.margin = margin(r = 1, l = 10, t = 10)),
                   p_scDEED+ 
                     theme(plot.margin = margin(r = 1, l = 1, t = 10)),
                   p_DynamicViz + 
                     theme(plot.margin = margin(r = 1, l = 10, t = 10)),nrow = 2)
ggsave(P, filename = './plots/Comparison/compare.png', width = 10/7*10, height = 6.8/7*10)

# Plot with spiral angle
p_spiral = ggplot(data = data.frame(x = Y[,1], y = Y[,2], t = t)) + geom_point(aes(x = x, y = y, color = t), size = 1.2) + xlab('tSNE1') + ylab('tSNE2')+ 
  theme_minimal()+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 16), 
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16)
  ) + labs(color = "Spiral\nAngle") + ggtitle('t-SNE Embedding with Spiral Angle')
ggsave(p_spiral, filename = './plots/Comparison/compare_spiralangle.png', width = 4.7/7*10, height = 3.4/7*10)

########################################
# Comparison with singularity score
########################################

########################################
# Mouse Brain

# EMBEDR for mouse brain data
library(RcppCNPy)
library(ggplot2)
library(reshape2)
perplexity_candidates = c(seq(5,150,5))
n = 3618
pvalues_mouse_brain = matrix(NA, nrow = n, ncol = length(perplexity_candidates))
for (i in c(1:length(perplexity_candidates))){
  pvalues_mouse_brain[,i] = npyLoad(paste0("./data/Comparison/embedr_mousebrain/",perplexity_candidates[i],".npy"))
}
pvalues_df <- as.data.frame(pvalues_mouse_brain)
colnames(pvalues_df) <- perplexity_candidates
pvalues_long <- melt(pvalues_df, variable.name = "Perplexity", value.name = "pvalue")
pvalues_long$CandidateType <- ifelse(as.character(pvalues_long$Perplexity) == as.character(145), "Minimal", "Other")
pvalues_long$Perplexity <- as.numeric(as.character(pvalues_long$Perplexity))
p_embedr = ggplot(pvalues_long, aes(x = Perplexity, y = pvalue, group = Perplexity, color = CandidateType)) +
  geom_boxplot(outlier.colour="darkgrey", outlier.shape=16,
               outlier.size=0.5, notch=TRUE)+
  scale_color_manual(values = c("Other" = "darkgrey", "Minimal" = "red"), guide = "none") +
  labs(title = "EMBEDR",
       x = "Perplexity",
       y = "EMBEDR p-value") +
  theme_minimal() + 
  scale_y_log10()+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16 ),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16),
    legend.position = "none"
  ) + scale_x_continuous(breaks = c(0, 50, 100, 150))

# scDEED for mouse brain data
load("./data/Comparison/scDEED_mousebrain.Rdata")
p_scDEED = ggplot() + geom_point(data = scDEED_result$num_dubious,aes(x = as.numeric(perplexity), y = number_dubious_cells)) + xlab("Perplexity") + ylab("Number of dubious cells") + labs(title = "scDEED")+
  theme_minimal()  + geom_point(aes(x = c(10,145), y = c(40,7)), color = 'red', size = 2)+
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 16 ), 
    plot.title = element_text(size = 20),  
    plot.subtitle = element_text(size = 20),  
    plot.caption = element_text(size = 14),  
    strip.text = element_text(size = 16),
    legend.position = "none"
  )

# DynamicViz for mouse brain data
perplexity_candidates = c(seq(5,150,5))
n = 3618
variancescore_mouse_brain = matrix(NA, nrow = n, ncol = length(perplexity_candidates))
for (i in c(1:length(perplexity_candidates))){
  variancescore_mouse_brain[,i] = read.csv(paste0("./data/Comparison/dynamicviz_mousebrain/dynamicviz_B20_tsne_perplexity",perplexity_candidates[i],"_variance.csv"), header = FALSE, sep = ",")[[1]]
}
var_df <- as.data.frame(variancescore_mouse_brain)
colnames(var_df) <- perplexity_candidates
var_long <- melt(var_df, variable.name = "Perplexity", value.name = "varianceScore")
var_long$Perplexity <- as.numeric(as.character(var_long$Perplexity))
calcStats <- function(x) {
  c(median = median(x),
    lower  = quantile(x, 0.025),
    upper  = quantile(x, 0.975))
}
summary_stats <- aggregate(varianceScore ~ Perplexity, data = var_long, FUN = calcStats)
stats <- data.frame(Perplexity = summary_stats$Perplexity, summary_stats$varianceScore)
colnames(stats) <- c("Perplexity", "median", "lower", "upper")
p_dynamicviz = ggplot(stats, aes(x = Perplexity)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgray", alpha = 0.5) +
  geom_line(aes(y = median), color = "black", size = 1, linetype = "dotted") +
  geom_point(aes(y = median), color = "black", size = 1.5) +
  labs(title = "DynamicViz",
       x = "Perplexity",
       y = "Variance Score") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 150)) + geom_point(data = data.frame(x = c(10), y = median(variancescore_mouse_brain[,2])), aes(x = x, y = y), color = 'red', size = 2)+
  theme(
    text = element_text(size = 16), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 16 ), 
    plot.title = element_text(size = 20), 
    plot.subtitle = element_text(size = 20), 
    plot.caption = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.position = "none"
  )

# Singularity score for mouse brain data
load("./data/Comparison/sscore_mousebrain.RData")
p_sscore = ggplot() + geom_point(aes(x = perplexity_candidates, y = mean_1_eigen_q)) + xlab("Perplexity") + ylab("Mean of top 5% singularity score") + labs(title = "Singularity Score")+theme_minimal()  + geom_point(aes(x = perplexity_candidates[19], y = mean_1_eigen_q[19]), color = "red", size = 2)+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16 ),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16),
    legend.position = "none"
  )
P_mouse_brain = egg::ggarrange(p_sscore, p_dynamicviz, p_scDEED, p_embedr, ncol = 4)
ggsave(P_mouse_brain, file = "./plots/Comparison/perplexity_selection_mouse_brain.png", width = 13.5/7*10, height = 3.6/7*10)

# plot embeddings
load("./data/Comparison/brain.RData")
label = as.factor(factor(ct.X1, labels = c('Astrocytes', 'Endothelial\nCells', 'Excitatory\nNeurons', 'Inhibitory\nNeurons', 'Microglia', 'Oligoden\n-drocytes')))
perplexity_vec = c(10, 30, 95, 145)
p_list = vector("list", 4)
i=1
filename = paste('./data/Comparison/small_atac_gene_activity_10geno_tsne_perplexity_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y1 = Y
i=2
filename = paste('./data/Comparison/small_atac_gene_activity_10geno_tsne_perplexity_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y2 = Y
i=3
filename = paste('./data/Comparison/small_atac_gene_activity_10geno_tsne_perplexity_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y3 = Y
i=4
filename = paste('./data/Comparison/small_atac_gene_activity_10geno_tsne_perplexity_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y4 = Y
P_embedding_mouse_brain = egg::ggarrange(ggplot() +
                                           geom_point(aes(x = Y1[,1], y = Y1[,2], color =label), size = 0.2, show.legend = TRUE)  +
                                           xlab('tSNE1') + ylab('tSNE2') +
                                           theme_minimal()+
                                           theme(
                                             text = element_text(size = 16),  # Increase font size for all text elements
                                             axis.title = element_text(size = 16), 
                                             axis.text = element_text(size = 14), 
                                             legend.title = element_text(size = 18), 
                                             legend.text = element_text(size = 16), 
                                             plot.title = element_text(size = 20),  # Increase font size for plot title
                                             plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                             plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                             strip.text = element_text(size = 20,hjust = 0)
                                           ) + ggtitle("Perplexity 10 (scDEED, DynamicViz)") +
                                           labs(color = "Cell Type") + 
                                           guides(color = guide_legend(override.aes = list(size = 2)))  + theme(legend.position = "none"),
                                         ggplot() +
                                           geom_point(aes(x = Y2[,1], y = Y2[,2], color =label), size = 0.2, show.legend = TRUE)  +
                                           xlab('tSNE1') + ylab('tSNE2') +
                                           theme_minimal()+
                                           theme(
                                             text = element_text(size = 16),  # Increase font size for all text elements
                                             axis.title = element_text(size = 16), 
                                             axis.text = element_text(size = 14), 
                                             legend.title = element_text(size = 18), 
                                             legend.text = element_text(size = 16), 
                                             plot.title = element_text(size = 20),  # Increase font size for plot title
                                             plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                             plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                             strip.text = element_text(size = 20,hjust = 0)
                                           ) +ggtitle("Perplexity 30 (Default)") +
                                           labs(color = "Cell Type") + 
                                           guides(color = guide_legend(override.aes = list(size = 2)))  + theme(axis.title.y = element_blank(),legend.position = "none"), 
                                         ggplot() +
                                           geom_point(aes(x = Y3[,1], y = Y3[,2], color =label), size = 0.2, show.legend = TRUE)  +
                                           xlab('tSNE1') + ylab('tSNE2') +
                                           theme_minimal()+
                                           theme(
                                             text = element_text(size = 16),  # Increase font size for all text elements
                                             axis.title = element_text(size = 16), 
                                             axis.text = element_text(size = 14), 
                                             legend.title = element_text(size = 18), 
                                             legend.text = element_text(size = 16), 
                                             plot.title = element_text(size = 20),  # Increase font size for plot title
                                             plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                             plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                             strip.text = element_text(size = 20,hjust = 0)
                                           ) +ggtitle("Perplexity 95 (Singularity Score)") +
                                           labs(color = "Cell Type") + 
                                           guides(color = guide_legend(override.aes = list(size = 2)))  + theme(axis.title.y = element_blank(),legend.position = "none"), 
                                         ggplot() +
                                           geom_point(aes(x = Y4[,1], y = Y4[,2], color =label), size = 0.2, show.legend = TRUE)  +
                                           xlab('tSNE1') + ylab('tSNE2') +
                                           theme_minimal()+
                                           theme(
                                             text = element_text(size = 16),  # Increase font size for all text elements
                                             axis.title = element_text(size = 16), 
                                             axis.text = element_text(size = 14), 
                                             legend.title = element_text(size = 18), 
                                             legend.text = element_text(size = 16), 
                                             plot.title = element_text(size = 20),  # Increase font size for plot title
                                             plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                             plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                             strip.text = element_text(size = 20,hjust = 0)
                                           ) +ggtitle("Perplexity 145 (scDEED, EMBEDR)") +
                                           labs(color = "Cell Type") + 
                                           guides(color = guide_legend(override.aes = list(size = 2)))  + theme(axis.title.y = element_blank()), ncol = 4)
ggsave(P_embedding_mouse_brain, file = "./plots/Comparison/embeddings_mouse_brain.png", width = 14.54/7*10, height = 3.6/7*10)

########################################
# Hayashi Embryo

# scDEED
load("./data/Comparison/scDEED_hayashi.RData")
perplexity = scDEED_result$num_dubious$perplexity
n_d_cells = scDEED_result$num_dubious$number_dubious_cells
p_scDEED = ggplot() + geom_point(data = scDEED_result$num_dubious,aes(x = as.numeric(perplexity), y = number_dubious_cells)) + xlab("Perplexity") + ylab("Number of dubious cells") + labs(title = "scDEED")+
  theme_minimal()  + geom_point(aes(x = as.numeric(c(3, perplexity[n_d_cells==1])), y = c(3, rep(1, sum(n_d_cells==1)))), color = 'red', size = 2)+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16 ),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16),
    legend.position = "none"
  )

# DynamicViz
perplexity_candidates = c(seq(2,40,1))
n = 421
variancescore_mouse_brain = matrix(NA, nrow = n, ncol = length(perplexity_candidates))
for (i in c(1:length(perplexity_candidates))){
  variancescore_mouse_brain[,i] = read.csv(paste0("./data/Comparison/dynamicviz_hayashi/dynamicviz_B100_tsne_perplexity",perplexity_candidates[i],"_variance.csv"), header = FALSE, sep = ",")[[1]]
}
var_df <- as.data.frame(variancescore_mouse_brain)
colnames(var_df) <- perplexity_candidates
var_long <- melt(var_df, variable.name = "Perplexity", value.name = "varianceScore")
var_long$Perplexity <- as.numeric(as.character(var_long$Perplexity))
calcStats <- function(x) {
  c(median = median(x),
    lower  = quantile(x, 0.025),
    upper  = quantile(x, 0.975))
}
summary_stats <- aggregate(varianceScore ~ Perplexity, data = var_long, FUN = calcStats)
stats <- data.frame(Perplexity = summary_stats$Perplexity, summary_stats$varianceScore)
colnames(stats) <- c("Perplexity", "median", "lower", "upper")
p_dynamicviz = ggplot(stats, aes(x = Perplexity)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgray", alpha = 0.5) +
  geom_line(aes(y = median), color = "black", size = 1, linetype = "dotted") +
  geom_point(aes(y = median), color = "black", size = 1.5) +
  labs(title = "DynamicViz",
       x = "Perplexity",
       y = "Variance Score") +
  theme_minimal() +
  geom_point(data = data.frame(x = c(20), y = median(variancescore_mouse_brain[,19])), aes(x = x, y = y), color = 'red', size = 2)+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16 ),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16),
    legend.position = "none"
  )

# Singularity score
load("./data/Comparison/sscore_hayashi.RData")
p_sscore = ggplot() + geom_point(aes(x = perplexity_candidates, y = mean_1_eigen_q)) + xlab("Perplexity") + ylab("Mean of top 5% singularity score") + labs(title = "Singularity Score")+theme_minimal()  + geom_point(aes(x = perplexity_candidates[24], y = mean_1_eigen_q[24]), color = "red", size = 2)+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16 ),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16),
    legend.position = "none"
  )
P_hayashi = egg::ggarrange(p_sscore, p_dynamicviz, p_scDEED, ncol = 3)
ggsave(P_hayashi, file = "./plots/Comparison/perplexity_selection_hayashi.png", width = 13.5/7*10/4*3, height = 3.6/7*10)

# plot embeddings
load("./data/Comparison/hayashi.RData")
perplexity_vec = c(3,20,25,30)
i=1
filename = paste('./data/Comparison/hayashi_tsne_perplexity_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y1 = Y
i=2
filename = paste('./data/Comparison/hayashi_tsne_perplexity_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y2 = Y
i=3
filename = paste('./data/Comparison/hayashi_tsne_perplexity_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y3 = Y
i=4
filename = paste('./data/Comparison/hayashi_tsne_perplexity_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y4 = Y
label = as.numeric(ct.X1)
label[label == 1] = '00 h'
label[label == 2] = '12 h'
label[label == 3] = '24 h'
label[label == 4] = '48 h'
label[label == 5] = '72 h'
P_embedding_hayashi = egg::ggarrange(ggplot() +
                                       geom_point(aes(x = Y1[,1], y = Y1[,2], color =label), size = 1, show.legend = TRUE)  +
                                       xlab('tSNE1') + ylab('tSNE2') +
                                       theme_minimal()+
                                       theme(
                                         text = element_text(size = 16),  # Increase font size for all text elements
                                         axis.title = element_text(size = 16), 
                                         axis.text = element_text(size = 14), 
                                         legend.title = element_text(size = 18), 
                                         legend.text = element_text(size = 16), 
                                         plot.title = element_text(size = 20),  # Increase font size for plot title
                                         plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                         plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                         strip.text = element_text(size = 20,hjust = 0)
                                       ) + ggtitle("Perplexity 3 (scDEED)") +
                                       labs(color = "Time\nPoint") + 
                                       guides(color = guide_legend(override.aes = list(size = 2)))  + theme(legend.position = "none"),
                                     ggplot() +
                                       geom_point(aes(x = Y2[,1], y = Y2[,2], color =label), size = 1, show.legend = TRUE)  +
                                       xlab('tSNE1') + ylab('tSNE2') +
                                       theme_minimal()+
                                       theme(
                                         text = element_text(size = 16),  # Increase font size for all text elements
                                         axis.title = element_text(size = 16), 
                                         axis.text = element_text(size = 14), 
                                         legend.title = element_text(size = 18), 
                                         legend.text = element_text(size = 16), 
                                         plot.title = element_text(size = 20),  # Increase font size for plot title
                                         plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                         plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                         strip.text = element_text(size = 20,hjust = 0)
                                       ) +ggtitle("Perplexity 20 (DynamicViz)") +
                                       labs(color = "Time\nPoint") + 
                                       guides(color = guide_legend(override.aes = list(size = 2)))  + theme(axis.title.y = element_blank(),legend.position = "none"), 
                                     ggplot() +
                                       geom_point(aes(x = Y3[,1], y = Y3[,2], color =label), size = 1, show.legend = TRUE)  +
                                       xlab('tSNE1') + ylab('tSNE2') +
                                       theme_minimal()+
                                       theme(
                                         text = element_text(size = 16),  # Increase font size for all text elements
                                         axis.title = element_text(size = 16), 
                                         axis.text = element_text(size = 14), 
                                         legend.title = element_text(size = 18), 
                                         legend.text = element_text(size = 16), 
                                         plot.title = element_text(size = 20),  # Increase font size for plot title
                                         plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                         plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                         strip.text = element_text(size = 20,hjust = 0)
                                       ) +ggtitle("Perplexity 25 (Singularity Score)") +
                                       labs(color = "Time\nPoint") + 
                                       guides(color = guide_legend(override.aes = list(size = 2)))  + theme(axis.title.y = element_blank(),legend.position = "none"), 
                                     ggplot() +
                                       geom_point(aes(x = Y4[,1], y = Y4[,2], color =label), size = 1, show.legend = TRUE)  +
                                       xlab('tSNE1') + ylab('tSNE2') +
                                       theme_minimal()+
                                       theme(
                                         text = element_text(size = 16),  # Increase font size for all text elements
                                         axis.title = element_text(size = 16), 
                                         axis.text = element_text(size = 14), 
                                         legend.title = element_text(size = 18), 
                                         legend.text = element_text(size = 16), 
                                         plot.title = element_text(size = 20),  # Increase font size for plot title
                                         plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                         plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                         strip.text = element_text(size = 20,hjust = 0)
                                       ) +ggtitle("Perplexity 30 (Default)") +
                                       labs(color = "Time\nPoint") + 
                                       guides(color = guide_legend(override.aes = list(size = 2)))  + theme(axis.title.y = element_blank()), ncol = 4)
ggsave(P_embedding_hayashi, file = "./plots/Comparison/embeddings_hayashi.png", width = 14.54/7*10, height = 3.6/7*10)

########################################
# Mammary

# EMBEDR
perplexity_candidates = c(seq(25,500,25))
n = 25806
pvalues_mouse_brain = matrix(NA, nrow = n, ncol = length(perplexity_candidates))
for (i in c(1:length(perplexity_candidates))){
  pvalues_mouse_brain[,i] = npyLoad(paste0("./data/Comparison/embedr_bachmammary/",perplexity_candidates[i],".npy"))
}
pvalues_df <- as.data.frame(pvalues_mouse_brain)
colnames(pvalues_df) <- perplexity_candidates
pvalues_long <- melt(pvalues_df, variable.name = "Perplexity", value.name = "pvalue")
pvalues_long$CandidateType <- ifelse(as.character(pvalues_long$Perplexity) == as.character(450), "Minimal", "Other")
pvalues_long$Perplexity <- as.numeric(as.character(pvalues_long$Perplexity))
p_embedr = ggplot(pvalues_long, aes(x = Perplexity, y = pvalue, group = Perplexity, color = CandidateType)) +
  geom_boxplot(outlier.colour="darkgrey", outlier.shape=16,
               outlier.size=0.5, notch=TRUE)+
  scale_color_manual(values = c("Other" = "darkgrey", "Minimal" = "red"), guide = "none") +
  labs(title = "EMBEDR",
       x = "Perplexity",
       y = "EMBEDR p-value") +
  theme_minimal() + 
  scale_y_log10()+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16 ),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16),
    legend.position = "none"
  )+
  scale_x_continuous(breaks = c(seq(0,500,100)))

# scDEED
scDEED_result = readRDS("./data/Comparison/scDEED_bachmammary.Rds")
p_scDEED = ggplot() + geom_point(data = scDEED_result$num_dubious,aes(x = as.numeric(perplexity), y = number_dubious_cells)) + xlab("Perplexity") + ylab("Number of dubious cells") + labs(title = "scDEED")+
  theme_minimal()  + geom_point(aes(x = c(100,500), y = c(303,129)), color = 'red', size = 2)+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16 ),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16),
    legend.position = "none"
  )

# Singularity score
load("./data/Comparison/sscore_bachmammary.RData")
mean_1_eigen_q = mean_1_eigen_q[-2]
p_sscore = ggplot() + geom_point(aes(x = perplexity_candidates, y = mean_1_eigen_q)) + xlab("Perplexity") + ylab("Mean of top 5% singularity score") + labs(title = "Singularity Score")+theme_minimal()  + geom_point(aes(x = perplexity_candidates[7], y = mean_1_eigen_q[7]), color = "red", size = 2)+
  theme(
    text = element_text(size = 16),  # Increase font size for all text elements
    axis.title = element_text(size = 16),  # Increase font size for axis titles
    axis.text = element_text(size = 14),  # Increase font size for axis text
    legend.title = element_text(size = 18),  # Increase font size for legend title
    legend.text = element_text(size = 16 ),  # Increase font size for legend text
    plot.title = element_text(size = 20),  # Increase font size for plot title
    plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
    plot.caption = element_text(size = 14),  # Increase font size for plot caption
    strip.text = element_text(size = 16),
    legend.position = "none"
  )
P_mouse_brain = egg::ggarrange(p_sscore, p_scDEED, p_embedr, ncol = 3)
ggsave(P_mouse_brain, file = "./plots/Comparison/perplexity_selection_bach_mammary.png", width = 13.5/7*10/4*3, height = 3.6/7*10)

# Plot embeddings
load("./data/Comparison/BachMammary_label.RData")
# label = as.factor(factor(bachmammary_label))
label = as.factor(factor(bachmammary_label,labels = c("G1", "G2", "L1", "L2", "NP1", "NP2", "PI1", "PI2")))
perplexity_vec = c(100, 175, 450, 500)
p_list = vector("list", 4)

i=1
filename = paste('./data/Comparison/BachMammaryData_tSNE_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y1 = Y
i=2
filename = paste('./data/Comparison/BachMammaryData_tSNE_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y2 = Y
i=3
filename = paste('./data/Comparison/BachMammaryData_tSNE_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y3 = Y
i=4
filename = paste('./data/Comparison/BachMammaryData_tSNE_',perplexity_vec[i],'.RData',sep='')
load(filename)
Y4 = Y
P_embedding_mouse_brain = egg::ggarrange(ggplot() +
                                           geom_point(aes(x = Y1[,1], y = Y1[,2], color =label), size = 0.2, show.legend = TRUE)  +
                                           xlab('tSNE1') + ylab('tSNE2') +
                                           theme_minimal()+
                                           theme(
                                             text = element_text(size = 16),  # Increase font size for all text elements
                                             axis.title = element_text(size = 16), 
                                             axis.text = element_text(size = 14), 
                                             legend.title = element_text(size = 18), 
                                             legend.text = element_text(size = 16), 
                                             plot.title = element_text(size = 20),  # Increase font size for plot title
                                             plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                             plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                             strip.text = element_text(size = 20,hjust = 0)
                                           ) + ggtitle("Perplexity 100 (scDEED)") +
                                           labs(color = "Cell Type") + 
                                           guides(color = guide_legend(override.aes = list(size = 2)))  + theme(legend.position = "none"),
                                         ggplot() +
                                           geom_point(aes(x = Y2[,1], y = Y2[,2], color =label), size = 0.2, show.legend = TRUE)  +
                                           xlab('tSNE1') + ylab('tSNE2') +
                                           theme_minimal()+
                                           theme(
                                             text = element_text(size = 16),  # Increase font size for all text elements
                                             axis.title = element_text(size = 16), 
                                             axis.text = element_text(size = 14), 
                                             legend.title = element_text(size = 18), 
                                             legend.text = element_text(size = 16), 
                                             plot.title = element_text(size = 20),  # Increase font size for plot title
                                             plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                             plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                             strip.text = element_text(size = 20,hjust = 0)
                                           ) +ggtitle("Perplexity 175 (Singularity Score)") +
                                           labs(color = "Cell Type") + 
                                           guides(color = guide_legend(override.aes = list(size = 2)))  + theme(axis.title.y = element_blank(),legend.position = "none"), 
                                         ggplot() +
                                           geom_point(aes(x = Y3[,1], y = Y3[,2], color =label), size = 0.2, show.legend = TRUE)  +
                                           xlab('tSNE1') + ylab('tSNE2') +
                                           theme_minimal()+
                                           theme(
                                             text = element_text(size = 16),  # Increase font size for all text elements
                                             axis.title = element_text(size = 16), 
                                             axis.text = element_text(size = 14), 
                                             legend.title = element_text(size = 18), 
                                             legend.text = element_text(size = 16), 
                                             plot.title = element_text(size = 20),  # Increase font size for plot title
                                             plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                             plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                             strip.text = element_text(size = 20,hjust = 0)
                                           ) +ggtitle("Perplexity 450 (EMBEDR)") +
                                           labs(color = "Cell Type") + 
                                           guides(color = guide_legend(override.aes = list(size = 2)))  + theme(axis.title.y = element_blank(),legend.position = "none"), 
                                         ggplot() +
                                           geom_point(aes(x = Y4[,1], y = Y4[,2], color =label), size = 0.2, show.legend = TRUE)  +
                                           xlab('tSNE1') + ylab('tSNE2') +
                                           theme_minimal()+
                                           theme(
                                             text = element_text(size = 16),  # Increase font size for all text elements
                                             axis.title = element_text(size = 16), 
                                             axis.text = element_text(size = 14), 
                                             legend.title = element_text(size = 18), 
                                             legend.text = element_text(size = 16), 
                                             plot.title = element_text(size = 20),  # Increase font size for plot title
                                             plot.subtitle = element_text(size = 20),  # Increase font size for plot subtitle
                                             plot.caption = element_text(size = 14),  # Increase font size for plot caption
                                             strip.text = element_text(size = 20,hjust = 0)
                                           ) +ggtitle("Perplexity 500 (scDEED)") +
                                           labs(color = "Cell Type") + 
                                           guides(color = guide_legend(override.aes = list(size = 2)))  + theme(axis.title.y = element_blank()), ncol = 4)
ggsave(P_embedding_mouse_brain, file = "./plots/Comparison/embeddings_bach_mammary.png", width = 14.54/7*10, height = 3.6/7*10)