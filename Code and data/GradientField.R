########################################
# May 13, 2025, Zhexuan Liu ############
# Gradient Field #######################
########################################

library(Rtsne)
library(ggplot2)
library(tsneMDBD)
library(Rfast)
library(uwot)

get_loss_pointwise_tsne=
  function (yy, Y, Ydist_sq, Mat) 
  {
    yy_Y_dist_sq = rowSums((Y - matrix(yy, nrow(Y), ncol(Y), 
                                       byrow = TRUE))^2)
    yydist_sq = rbind(matrix(yy_Y_dist_sq, nrow = 1), as.matrix(Ydist_sq))
    yydist_sq = cbind(matrix(c(0, yy_Y_dist_sq), ncol = 1), yydist_sq)
    QMat = (1 + yydist_sq)^(-1)
    QMat = QMat/(sum(QMat) - dim(yydist_sq)[1])
    eps = 2^(-52)
    PMat = Mat + eps
    QMat = QMat + eps
    # return(sum(PMat * log(PMat/QMat)))
    ind = 1
    return(- 2 * sum(PMat[ind,] * log(QMat[ind,])) + log(sum(colsums(QMat))))
  }
gradient_compute_pointwise = function(yy, Y, P){
  Y = c(t(rbind(yy,Y)))
  n = dim(P)[1]
  n2 = dim(P)[2]
  n_Y = length(Y)
  if ((n!=n2) | (2*n!=n_Y)){
    stop('Input dimensions do not match!')
  }
  
  YMat = matrix(Y, 2, n) # make Y 2*n matrix
  YDiffDistSq = as.matrix(dist(t(YMat))**2)
  YDiffDistSq_double = rbind(YDiffDistSq[1,],YDiffDistSq[1,])
  P_double = rbind(P[1,],P[1,])
  deno = sum(1 / (1 + YDiffDistSq)) - n
  
  yy_double = matrix(rep(yy,n),2,n)
  yy_Y_diff = yy_double - YMat
  I1 = 4 * P_double * yy_Y_diff / (1 + YDiffDistSq_double)
  I2 = (-1) * 4 * yy_Y_diff / deno / ((1 + YDiffDistSq_double) ^ 2) * sum(P) # sum(P) = 1, but I just want to put it here for inner peace for now...
  
  G = rowSums(I1 + I2)
  
  return(G)
}

load("./data/2GMM/2d_2GMM_perplexity50.RData")
perplexity = 50
mean_1 = colMeans(X[label =='1',])
mean_2 = colMeans(X[label =='2',])

angle = 38/180
Y[,2] = Y[,2] +1.2
Y[,1] = Y[,1] 
Y[,1] = Y[,1]*cos(angle)
Y[,2] = Y[,2]*sin(angle)

perplexity = 50
mean_1 = colMeans(X[label =='1',])
mean_2 = colMeans(X[label =='2',])
dens = 21
# Define a scalar function
X_new = rbind(0.525* (mean_2 - mean_1) + mean_1, X)
tsne.result = Rtsne(X_new, perplexity = perplexity, theta = 0, 
                    max_iter = 1, check_duplicates = FALSE)
PMat = tsne.result$P

colnames(Y) = c('x', 'y')
r = 4
xr = c(-r,
       r) # margin for x axis
yr = c(-r,
       r) # margin for y axis

xr_plot = xr
xr_plot[1] = xr_plot[1]
xr_plot[2] = xr_plot[2]
yr_plot = yr
yr_plot[1] = yr_plot[1]
yr_plot[2] = yr_plot[2] 

### create plot grid
xgrid = c(0:dens) / dens * (xr[2] - xr[1]) + xr[1]
ygrid = c(0:dens) / dens * (yr[2] - yr[1]) + yr[1]
plotdf = data.frame(x = c(rep(xgrid, each = length(ygrid))),
                    y = c(rep(ygrid, length(xgrid))),
                    z = c(rep(0, length(xgrid) * length(ygrid))),
                    dx = c(rep(0, length(xgrid) * length(ygrid))),
                    dy = c(rep(0, length(xgrid) * length(ygrid))),
                    dxa = c(rep(0, length(xgrid) * length(ygrid))),
                    dya = c(rep(0, length(xgrid) * length(ygrid))))

Ydist_sq = as.matrix(dist(Y)**2)
theta = (colmeans(Y[label =='1',]) - colmeans(Y[label =='2',]))/2
for (i in c(1:dim(plotdf)[1])){
  y1 = plotdf[i,1]
  y2 = plotdf[i,2]
  plotdf[i,3] = get_loss_pointwise_tsne(c(y1,y2), Y, Ydist_sq, PMat)
  plotdf[i,c(4,5)] = gradient_compute_pointwise(c(y1,y2), Y, PMat)#gradient_compute(rbind(c(y1,y2), Y), PMat)[1:2]
  y_parallel = (sum(theta * c(y1,y2))/sum(theta^2)) * theta
  y_prop = c(y1,y2) - y_parallel
  plotdf[i,c(6,7)] = (y_parallel - y_prop) / sum(theta**2)
}


P2 = ggplot(data = plotdf)  +
  geom_contour(aes(x = x, y = y, z = as.numeric(z), colour = after_stat(level)), alpha = 0.5, bins = 7)+
  scale_color_viridis_c(direction = 1) +
  geom_segment(aes(x=x,y=y,xend = x - dx, yend = y - dy), 
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "orange") +
  guides(color = guide_colorbar(title = "t-SNE Loss"))+
  theme(legend.key.size = unit(4, "lines"))+
  xlab('tSNE1') +
  ylab('tSNE2') +
  theme_minimal()+
  theme(
    text = element_text(family = "Arial",size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    plot.title = element_text(size = 19),
    plot.subtitle = element_text(size = 19),
    plot.caption = element_text(size = 15),
    strip.text = element_text(size = 15)
  ) + ggtitle('Gradient Field of Real Loss')
P3 = ggplot(data = plotdf)  +
  geom_contour(aes(x = x, y = y, z = as.numeric(z), colour = after_stat(level)), alpha = 0.5, bins = 7)+
  scale_color_viridis_c(direction = 1) +
  geom_segment(aes(x=x,y=y,xend = x + dxa, yend = y + dya), 
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "orange") +
  guides(color = guide_colorbar(title = "t-SNE Loss"))+
  theme(legend.key.size = unit(4, "lines"))+
  xlab('tSNE1') +
  ylab('tSNE2') +
  theme_minimal()+
  theme(
    text = element_text(family = "Arial",size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    plot.title = element_text(size = 19),
    plot.subtitle = element_text(size = 19),
    plot.caption = element_text(size = 15),
    strip.text = element_text(size = 15)
  ) + ggtitle('Theoretical Gradient Field')

# ggsave(P1, file = '2_10_1.png', width = 6.5/7*10, height = 4.5/7*10)
ggsave(P2 + theme(legend.position = 'none'), file = './plots/GradientField/Gradient Field of Real Loss.png', width = 5, height = 4.1)
ggsave(P3 + theme(legend.position = 'none'), file = './plots/GradientField/Theoretical Gradient Field.png', width = 5, height = 4.1)

pdf("./plots/GradientField/figure8a.pdf", width = 52.5*3 / 25.4, height = 40*3 / 25.4)
P2+ theme(legend.position = 'none')
dev.off()

pdf("./plots/GradientField/figure8b.pdf", width = 52.5*3 / 25.4, height = 40*3 / 25.4)
P3+ theme(legend.position = 'none')
dev.off()
