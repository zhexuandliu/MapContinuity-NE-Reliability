library(Rfast)
library(viridis)
get_loss_pointwise = function(yy, ind, Mat, method = 'tsne', a=1){
  if (method == 'umap'){
    yyMat = t(yy)
    newpoint =  yyMat[,ind]
    yydistvec = yyMat - newpoint
    yyDiffDistSq = colsums(yydistvec ** 2)[-ind]
    I1 = sum(log(1 + a * yyDiffDistSq))
    vvec = Mat[-ind, ind]
    I2 = -Crossprod(as.matrix(1 - vvec), as.matrix(log(yyDiffDistSq)))[1, 1]
    return(I1 + I2)
  }else if (method == 'largevis'){
    nn = Mat[1,]$nn$euclidean$idx[ind,]
    Mat = Mat[2,]$similarity_graph
    yyMat = t(yy)
    newpoint =  yyMat[,ind]
    yydistvec = yyMat - newpoint
    yyDiffDistSq = colsums(yydistvec ** 2)[-ind]
    w = 1 / (1+yyDiffDistSq)
    v = c(Mat[-ind,ind])
    I1 = sum(v * log(w))
    I2 = 7 * sum(log(1-w[v==0]))
    return(I1+I2)
  }
  else{
    QMat = as.matrix(1+Dist(yy, square = TRUE))
    # QMat = QMat / (sum(colsums(QMat)) - dim(yy)[1])
    loss = 2 * sum(Mat[ind,] * log(QMat[ind,])) + log(sum(colsums(1/QMat))-dim(QMat)[1])
    return(loss)
  }
}

library(uwot)
getP = function(X, perplexity, method = 'tsne'){
  if (method == 'umap'){
    P = as.matrix(similarity_graph(X, n_neighbors = perplexity))
    return(P)
  }else if(method == 'largevis'){
    P = as.matrix(similarity_graph(X, n_neighbors = perplexity, method = 'largevis', ret_extra = c('nn')))
    return(P)
  }else{
    Y = Rtsne(X, perplexity = perplexity, theta = 0, max_iter = 1)
    return(Y$P)
  }
}

plot_interpolate = function(X, perplexity, Y, mean_1, mean_2, label, inter_dens, ind_vec, i, size_all = 0.5, size_min = 3, method = 'tsne', dens = 10, thr){

  colnames(Y) = c('x','y')
  xr = c(min(Y[,1]),
         max(Y[,1])) # margin for x axis
  yr = c(min(Y[,2]),
         max(Y[,2]) + (max(Y[,2]) - min(Y[,2]))/10) # margin for y axis

  xr_plot = xr
  xr_plot[1] = xr_plot[1]-0.01
  xr_plot[2] = xr_plot[2]+0.01
  yr_plot = yr
  yr_plot[1] = yr_plot[1]-0.01
  yr_plot[2] = yr_plot[2]+0.01

  ### create plot grid
  xgrid = c(0:dens) / dens * (xr[2] - xr[1]) + xr[1]
  ygrid = c(0:dens) / dens * (yr[2] - yr[1]) + yr[1]
  plotdf = data.frame(x = c(rep(xgrid, each = length(ygrid))),
                      y = c(rep(ygrid, length(xgrid))),
                      z = c(rep(0, length(xgrid) * length(ygrid))))

  ind = ind_vec[i]
  X_int = ind / (inter_dens + 1) * (mean_2 - mean_1) + mean_1 # interpolate from mean_1 to mean_2
  X_new = rbind(X, X_int)
  labelnew = c(label, 'new')
  P_new = getP(X_new, perplexity, method)

  library(parallel)
  no_cores = detectCores() - 1

  plotdf$z = unlist(mclapply(1:dim(plotdf)[1],
                             function(i)
                               return(get_loss_pointwise(rbind(Y, c(plotdf$x[i], plotdf$y[i])),
                                                         dim(Y)[1]+1,
                                                         P_new,
                                                         method)),
                             mc.cores = no_cores))

  min_ind = which.min(plotdf$z)
  plotdf$z[plotdf$z>thr] = thr
  #, breaks = c(seq(1000,9000,1000))
  p = ggplot(data = plotdf)  +
    geom_contour(aes(x = x, y = y, z = as.numeric(z), colour = after_stat(level)), bins = 8, size = 1)
  p = p+
    geom_point(data = data.frame(Y[label == unique(label)[1],],label[label == unique(label)[1]]), aes(x = x, y = y), color = '#009ADE', size = size_all, alpha = 0.25)+
    geom_point(data = data.frame(Y[label == unique(label)[2],],label[label == unique(label)[2]]), aes(x = x, y = y), color = '#FF1F5B', size = size_all, alpha = 0.25)+
    #scale_color_manual(values = c("red", "blue")) +
    geom_point(data = plotdf[min_ind, ], aes(x = x, y = y), color = 'orange', size = size_min, shape = 17) +
    xlab(NULL) +
    #theme_bw() +
    # theme(panel.grid = element_blank()) +
    theme(panel.border = element_blank()) +
    scale_color_viridis_c(direction = 1) +
    xlim(xr_plot) +
    ylim(yr_plot) +
    xlab('tSNE1') +
    ylab('tSNE2') +
    guides(color = guide_colorbar(title = "LOO Loss"))+
    theme(legend.key.size = unit(2, "lines"))+
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

  ########## draw bar ##########
  rectangle_data = data.frame(xmin =  -(max(plotdf$x) - min(plotdf$x)) / 4,
                              xmax =  +(max(plotdf$x) - min(plotdf$x)) / 4,
                              ymin = 0,
                              ymax = (max(plotdf$x) - min(plotdf$x)) / 2 / 15)
  rectangle_data2 = data.frame(xmin = (-15 + 29/(inter_dens-1)*(ind-1)) * (max(plotdf$x) - min(plotdf$x)) / 60,
                               xmax = (-14 + 29/(inter_dens-1)*(ind-1)) * (max(plotdf$x) - min(plotdf$x)) / 60,
                               ymin =   (max(plotdf$x) - min(plotdf$x)) / 2 / 15 / 4,
                               ymax =   (max(plotdf$x) - min(plotdf$x)) / 2 / 15 / 4 * 3)
  custom_labels = data.frame(
    x = c(0, inter_dens + 1),
    y = c(-0.5, -0.5),
    label = c("start point", "end point")
  )
  bar2 = ggplot() +
    geom_rect(data = rectangle_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", color = "black") +
    geom_rect(data = rectangle_data2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "red", color = "black") +
    theme_bw() +
    theme(plot.margin = margin(0, 0, 0, 0),
          legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 10, hjust = 0.5)) +
    coord_fixed(ratio = 1) +
    xlab(NULL)
  ############################################################
  p = p +
    annotation_custom(
      ggplotGrob(bar2),
      xmin = (max(plotdf$x) + min(plotdf$x)) / 2 - (max(plotdf$x) - min(plotdf$x)) / 4,
      xmax = (max(plotdf$x) + min(plotdf$x)) / 2 + (max(plotdf$x) - min(plotdf$x)) / 4,
      ymin = max(plotdf$y)
    )
  # if (i %in% c(1,2,3,4,5,6)){
  #   p = p +
  #     theme(axis.text.x = element_blank(),
  #           axis.line.x = element_blank(),
  #           axis.ticks.x = element_blank())
  # }
  #
  # if (i %in% c(2,3,5,6,8,9)){
  #   p = p +
  #     theme(axis.text.y = element_blank(),
  #           axis.line.y = element_blank(),
  #           axis.ticks.y = element_blank())
  # }
  print(c(min(plotdf$z),max(plotdf$z)))
  return(p)
}



plot_interpolate_9_plots = function(X, perplexity, Y, mean_1, mean_2, label, inter_dens, ind_vec, size_all = 0.5, size_min = 3, method = 'tsne', dens = 10, thr){
  library(patchwork)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(Rfast)
  P = egg::ggarrange(plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                                 ind_vec, 1, size_all = size_all, size_min = size_min, method, dens = dens, thr) +
                  theme(plot.margin = margin(r = 1),
                        legend.position = 'none'),
                plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                                 ind_vec, 2, size_all = size_all, size_min = size_min, method, dens = dens, thr) +
                  theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank(),
                        plot.margin = margin(r = 1, l = 1),
                        legend.position = 'none'),
                plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                                 ind_vec, 3, size_all = size_all, size_min = size_min, method, dens = dens, thr) +
                  theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank(),
                        plot.margin = margin(r = 1, l = 1),
                        legend.position = 'none'),
                plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                                 ind_vec, 4, size_all = size_all, size_min = size_min, method, dens = dens, thr)+
                  theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank(),
                        plot.margin = margin(r = 1, l = 1) ),
                # plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                #                  ind_vec, 5, size_all = size_all, size_min = size_min, method, dens = dens),
                # plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                #                  ind_vec, 6, size_all = size_all, size_min = size_min, method, dens = dens),
                # plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                #                  ind_vec, 7, size_all = size_all, size_min = size_min, method, dens = dens),
                # plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                #                  ind_vec, 8, size_all = size_all, size_min = size_min, method, dens = dens),
                # plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                #                  ind_vec, 9, size_all = size_all, size_min = size_min, method, dens = dens),
                nrow = 1)
  return(P)
}

plot_interpolate_gif_plots = function(X, perplexity, Y, mean_1, mean_2, label, inter_dens, ind_vec, size_all = 0.5, size_min = 3, method = 'tsne', dens = 10){
  library(patchwork)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(Rfast)
  for (i in ind_vec){
    filename = paste('landscape-ifnb-',i,'.pdf',sep = '')
    ggsave(filename,plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                                     ind_vec, i, size_all = size_all, size_min = size_min, method, dens = dens) +
             theme(legend.position = "none"), width = 8, height = 6)
  }
}



plot_interpolate_2dori = function(X, perplexity, Y, mean_1, mean_2, label, inter_dens, ind_vec, i, size_all = 0.5, size_min = 3, method = 'tsne', dens = 10){

  colnames(Y) = c('x','y')
  xr = c(min(Y[,1]),
         max(Y[,1])) # margin for x axis
  yr = c(min(Y[,2]),
         max(Y[,2]) + (max(Y[,2]) - min(Y[,2]))/10) # margin for y axis

  xr_plot = xr
  xr_plot[1] = xr_plot[1]-0.01
  xr_plot[2] = xr_plot[2]+0.01
  yr_plot = yr
  yr_plot[1] = yr_plot[1]-0.01
  yr_plot[2] = yr_plot[2]+0.01

  ### create plot grid
  xgrid = c(0:dens) / dens * (xr[2] - xr[1]) + xr[1]
  ygrid = c(0:dens) / dens * (yr[2] - yr[1]) + yr[1]
  plotdf = data.frame(x = c(rep(xgrid, each = length(ygrid))),
                      y = c(rep(ygrid, length(xgrid))),
                      z = c(rep(0, length(xgrid) * length(ygrid))))

  ind = ind_vec[i]
  X_int = ind / (inter_dens + 1) * (mean_2 - mean_1) + mean_1 # interpolate from mean_1 to mean_2
  X_new = rbind(X, X_int)
  labelnew = c(label, 'new')

  min_ind = which.min(plotdf$z)
  #, breaks = c(seq(1000,9000,1000))
  p = ggplot()
  p = p+
    geom_point(data = data.frame(Y,label), aes(x = x, y = y, color = label), size = size_all, alpha = 0.25)+
    scale_color_manual(values = c('#FF1F5B','#009ADE')) +
    geom_point(data = data.frame(x = X_int[1], y = X_int[2]), aes(x = x, y = y), color = 'orange', size = size_min) +
    xlab(NULL) +
    xlim(xr_plot) +
    ylim(yr_plot) +
    xlab('X1') +
    ylab('X2')+
    theme_minimal()+
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

  ########## draw bar ##########
  rectangle_data = data.frame(xmin =  -(max(plotdf$x) - min(plotdf$x)) / 4,
                              xmax =  +(max(plotdf$x) - min(plotdf$x)) / 4,
                              ymin = 0,
                              ymax = (max(plotdf$x) - min(plotdf$x)) / 2 / 15)
  rectangle_data2 = data.frame(xmin = (-15 + 29/(inter_dens-1)*(ind-1)) * (max(plotdf$x) - min(plotdf$x)) / 60,
                               xmax = (-14 + 29/(inter_dens-1)*(ind-1)) * (max(plotdf$x) - min(plotdf$x)) / 60,
                               ymin =   (max(plotdf$x) - min(plotdf$x)) / 2 / 15 / 4,
                               ymax =   (max(plotdf$x) - min(plotdf$x)) / 2 / 15 / 4 * 3)
  custom_labels = data.frame(
    x = c(0, inter_dens + 1),
    y = c(-0.5, -0.5),
    label = c("start point", "end point")
  )
  bar2 = ggplot() +
    geom_rect(data = rectangle_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", color = "black") +
    geom_rect(data = rectangle_data2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "red", color = "black") +
    theme_bw() +
    theme(plot.margin = margin(0, 0, 0, 0),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 10, hjust = 0.5)) +
    coord_fixed(ratio = 1) +
    xlab(NULL)
  ############################################################
  p = p +
    annotation_custom(ggplotGrob(bar2),
                      xmin = (max(plotdf$x) + min(plotdf$x)) / 2 - (max(plotdf$x) - min(plotdf$x)) / 4,
                      xmax = (max(plotdf$x) + min(plotdf$x)) / 2 + (max(plotdf$x) - min(plotdf$x)) / 4,
                      ymin = max(plotdf$y))
  return(p)
}

plot_interpolate_9_plots_2dori = function(X, perplexity, Y, mean_1, mean_2, label, inter_dens, ind_vec, size_all = 0.5, size_min = 3, method = 'tsne', dens = 10){
  library(patchwork)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(Rfast)
  P = egg::ggarrange(plot_interpolate_2dori(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                                 ind_vec, 1, size_all = size_all, size_min = size_min, method, dens = dens) +
                       theme(plot.margin = margin(r = 1),
                             legend.position = 'none'),
                plot_interpolate_2dori(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                                 ind_vec, 2, size_all = size_all, size_min = size_min, method, dens = dens) +
                  theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank(),
                        plot.margin = margin(r = 1, l = 1),
                        legend.position = 'none'),
                plot_interpolate_2dori(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                                 ind_vec, 3, size_all = size_all, size_min = size_min, method, dens = dens) +
                  theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank(),
                        plot.margin = margin(r = 1, l = 1),
                        legend.position = 'none'),
                plot_interpolate_2dori(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                                 ind_vec, 4, size_all = size_all, size_min = size_min, method, dens = dens) +
                  theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank(),
                        plot.margin = margin(r = 1, l = 1) ),
                # plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                #                  ind_vec, 5, size_all = size_all, size_min = size_min, method, dens = dens),
                # plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                #                  ind_vec, 6, size_all = size_all, size_min = size_min, method, dens = dens),
                # plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                #                  ind_vec, 7, size_all = size_all, size_min = size_min, method, dens = dens),
                # plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                #                  ind_vec, 8, size_all = size_all, size_min = size_min, method, dens = dens),
                # plot_interpolate(X, perplexity, Y, mean_1, mean_2, label, inter_dens,
                #                  ind_vec, 9, size_all = size_all, size_min = size_min, method, dens = dens),
                nrow = 1)
  return(P)
}



speed_eigenscore_interpolate = function(X, perplexity, Y, mean_1, mean_2, label, inter_dens, ind_vec, i, size_all = 0.5, size_min = 3, method = 'tsne', dens = 10){

  colnames(Y) = c('x','y')
  xr = c(min(Y[,1]) - (max(Y[,1])-min(Y[,1]))/4,
         max(Y[,1]) + (max(Y[,1])-min(Y[,1]))/4) # margin for x axis
  yr = c(min(Y[,2]) - (max(Y[,2])-min(Y[,2]))/4,
         max(Y[,2]) + (max(Y[,2]) - min(Y[,2]))/4) # margin for y axis

  ### create plot grid
  xgrid = c(0:dens) / dens * (xr[2] - xr[1]) + xr[1]
  ygrid = c(0:dens) / dens * (yr[2] - yr[1]) + yr[1]
  plotdf = data.frame(x = c(rep(xgrid, each = length(ygrid))),
                      y = c(rep(ygrid, length(xgrid))),
                      z = c(rep(0, length(xgrid) * length(ygrid))))

  ind = ind_vec[i]
  X_int = ind / (inter_dens + 1) * (mean_2 - mean_1) + mean_1 # interpolate from mean_1 to mean_2
  X_new = rbind(X, X_int)
  labelnew = c(label, 'new')

  Y_new = Rtsne(X_new, perplexity = perplexity, theta = 0, max_iter = 10000)
  eigen_score_new = eigen_score(Y_new$Y, Y_new$P)[dim(X_new)[1]]
  P_new = Y_new$P

  library(parallel)
  no_cores = detectCores() - 1
  plotdf$z = unlist(mclapply(1:dim(plotdf)[1],
                             function(i)
                               return(get_loss_pointwise(rbind(Y, c(plotdf$x[i], plotdf$y[i])),
                                                         dim(Y)[1]+1,
                                                         P_new,
                                                         method)),
                             mc.cores = no_cores))

  min_ind = which.min(plotdf$z)


  #####################################
  return(list(coord = c(plotdf$x[min_ind],plotdf$y[min_ind]), score = eigen_score_new))
}

coord_min_interpolate = function(X, perplexity, Y, mean_1, mean_2, label, inter_dens, ind_vec, i, method = 'tsne', dens = 10){

  colnames(Y) = c('x','y')
  xr = c(min(Y[,1]) - (max(Y[,1])-min(Y[,1]))/10,
         max(Y[,1]) + (max(Y[,1])-min(Y[,1]))/10) # margin for x axis
  yr = c(min(Y[,2]) - (max(Y[,2])-min(Y[,2]))/10,
         max(Y[,2]) + (max(Y[,2]) - min(Y[,2]))/10) # margin for y axis

  ### create plot grid
  xgrid = c(0:dens) / dens * (xr[2] - xr[1]) + xr[1]
  ygrid = c(0:dens) / dens * (yr[2] - yr[1]) + yr[1]
  plotdf = data.frame(x = c(rep(xgrid, each = length(ygrid))),
                      y = c(rep(ygrid, length(xgrid))),
                      z = c(rep(0, length(xgrid) * length(ygrid))))

  ind = ind_vec[i]
  X_int =   mean_1 * (1-ind/inter_dens) + mean_2 * ind/inter_dens # interpolate from mean_1 to mean_2
  X_new = rbind(X, X_int)
  labelnew = c(label, 'new')

  Y_new = Rtsne(X_new, perplexity = perplexity, theta = 0, max_iter = 1)
  P_new = Y_new$P

  library(parallel)
  no_cores = detectCores() - 1
  plotdf$z = unlist(mclapply(1:dim(plotdf)[1],
                             function(i)
                               return(get_loss_pointwise(rbind(Y, c(plotdf$x[i], plotdf$y[i])),
                                                         dim(Y)[1]+1,
                                                         P_new,
                                                         method)),
                             mc.cores = no_cores))

  min_ind = which.min(plotdf$z)


  #####################################
  return(c(plotdf$x[min_ind],plotdf$y[min_ind]))
}

## no interpolation
plot_no_interpolate = function(X, perplexity, Y, mean_1, mean_2, label, inter_dens, ind_vec, i, size_all = 0.5, size_min = 3, method = 'tsne', dens = 10, thr){
  
  colnames(Y) = c('x','y')
  xr = c(min(Y[,1]),
         max(Y[,1])) # margin for x axis
  yr = c(min(Y[,2]),
         max(Y[,2]) + (max(Y[,2]) - min(Y[,2]))/10) # margin for y axis
  
  xr_plot = xr
  xr_plot[1] = xr_plot[1]-0.01
  xr_plot[2] = xr_plot[2]+0.01
  yr_plot = yr
  yr_plot[1] = yr_plot[1]-0.01
  yr_plot[2] = yr_plot[2]+0.01
  
  ### create plot grid
  xgrid = c(0:dens) / dens * (xr[2] - xr[1]) + xr[1]
  ygrid = c(0:dens) / dens * (yr[2] - yr[1]) + yr[1]
  plotdf = data.frame(x = c(rep(xgrid, each = length(ygrid))),
                      y = c(rep(ygrid, length(xgrid))),
                      z = c(rep(0, length(xgrid) * length(ygrid))))
  
  ind = ind_vec[i]
  X_int = ind / (inter_dens + 1) * (mean_2 - mean_1) + mean_1 # interpolate from mean_1 to mean_2
  X_new = rbind(X, X_int)
  labelnew = c(label, 'new')
  
  
  P_new = getP(X_new, perplexity, method)
  
  library(parallel)
  no_cores = detectCores() - 1
  
  plotdf$z = unlist(mclapply(1:dim(plotdf)[1],
                             function(i)
                               return(get_loss_pointwise(rbind(Y, c(plotdf$x[i], plotdf$y[i])),
                                                         dim(Y)[1]+1,
                                                         P_new,
                                                         method)),
                             mc.cores = no_cores))
  
  min_ind = which.min(plotdf$z)
  plotdf$z[plotdf$z>thr] = thr
  #, breaks = c(seq(1000,9000,1000))
  p = ggplot(data = plotdf)  +
    geom_contour(aes(x = x, y = y, z = as.numeric(z), colour = after_stat(level)), bins = 8, size = 1)
  p = p+
    geom_point(data = data.frame(x = Y[,1], y = Y[,2]), aes(x = x, y = y), color = 'black', size = size_all, alpha = 0.25)+
    geom_point(data = plotdf[min_ind, ], aes(x = x, y = y), color = 'orange', size = size_min, shape = 17) +
    xlab(NULL) +
    #theme_bw() +
    # theme(panel.grid = element_blank()) +
    theme(panel.border = element_blank()) +
    scale_color_viridis_c(direction = 1) +
    xlim(xr_plot) +
    ylim(yr_plot) +
    xlab('tSNE1') +
    ylab('tSNE2') +
    guides(color = guide_colorbar(title = "LOO Loss"))+
    theme(legend.key.size = unit(2, "lines"))+
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
  
  ########## draw bar ##########
  rectangle_data = data.frame(xmin =  -(max(plotdf$x) - min(plotdf$x)) / 4,
                              xmax =  +(max(plotdf$x) - min(plotdf$x)) / 4,
                              ymin = 0,
                              ymax = (max(plotdf$x) - min(plotdf$x)) / 2 / 15)
  rectangle_data2 = data.frame(xmin = (-15 + 29/(inter_dens-1)*(ind-1)) * (max(plotdf$x) - min(plotdf$x)) / 60,
                               xmax = (-14 + 29/(inter_dens-1)*(ind-1)) * (max(plotdf$x) - min(plotdf$x)) / 60,
                               ymin =   (max(plotdf$x) - min(plotdf$x)) / 2 / 15 / 4,
                               ymax =   (max(plotdf$x) - min(plotdf$x)) / 2 / 15 / 4 * 3)
  custom_labels = data.frame(
    x = c(0, inter_dens + 1),
    y = c(-0.5, -0.5),
    label = c("start point", "end point")
  )
  bar2 = ggplot() +
    geom_rect(data = rectangle_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", color = "black") +
    geom_rect(data = rectangle_data2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "red", color = "black") +
    theme_bw() +
    theme(plot.margin = margin(0, 0, 0, 0),
          legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 10, hjust = 0.5)) +
    coord_fixed(ratio = 1) +
    xlab(NULL)
  ############################################################
  p = p 
  # if (i %in% c(1,2,3,4,5,6)){
  #   p = p +
  #     theme(axis.text.x = element_blank(),
  #           axis.line.x = element_blank(),
  #           axis.ticks.x = element_blank())
  # }
  #
  # if (i %in% c(2,3,5,6,8,9)){
  #   p = p +
  #     theme(axis.text.y = element_blank(),
  #           axis.line.y = element_blank(),
  #           axis.ticks.y = element_blank())
  # }
  print(c(min(plotdf$z),max(plotdf$z)))
  return(p)
}