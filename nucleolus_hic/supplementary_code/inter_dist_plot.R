#################################################################
# NAIR 区的region plot inter-dist  但是使用的是已经算好的RData
#################################################################
# 作为补充材料，将其画到1张fig 中

rm(list=ls())

binsize = 100e3

fig.facet <- layout(matrix(c(1:24),nrow = 6,byrow = FALSE),width = rep(10,24),heights = rep(10,24))
layout.show(fig.facet)


for(chrom_index in c(1:23)){
  print(chrom_index)
  load(sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela_inter_dist/hela_inter_dist_chr_%d_%d.RData",chrom_index,binsize))
  
  ylim = c(0,7)
  xlim = c(5,8.5)
  
  par(mar=c(2,2,2,1))
  plot(x=log10(g_hic_dist_table$start),y=log10(g_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#1533AD",frame.plot=F,yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  axis(side=2,at=c(0:7))
  
  par(new=T)
  plot(x=log10(g_hic_ActD_dist_table$start),y=log10(g_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#6F81D6",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  par(new=T)
  plot(x=log10(n_hic_dist_table$start),y=log10(n_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FF8B00",frame.plot=F,yaxt="n",xlab="",xaxt="n",ylab="",xlim=xlim,ylim=ylim)
  
  par(new=T)
  plot(x=log10(n_hic_ActD_dist_table$start),y=log10(n_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FFC640",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  title(main = sprintf("chr%d NAIR",chrom_index))
  
}


###################################################################
# plot Hela-g-hic & Hela-g-hic-ActD interction-distance figure 
###################################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/make_dist_table.RData")
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")

binsize = 100e3

fig.facet <- layout(matrix(c(1:24),nrow = 6,byrow = FALSE),width = rep(10,24),heights = rep(10,24))
layout.show(fig.facet)

for(chrom_index in c(1:23)){
  print(chrom_index)
  
  load(file=sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_all_matrix_chr_%d_%d.RData",chrom_index,binsize))
  
  # conver as.matrix
  g_hic_matrix.raw = as.matrix(g_hic_matrix)
  g_hic_ActD_matrix.raw = as.matrix(g_hic_ActD_matrix)
  
  # make matrix with the same sequencing depth
  g_hic_matrix.fix = g_hic_matrix.raw / hic_count.hela_g_hic.total * 300e6
  g_hic_ActD_matrix.fix = g_hic_ActD_matrix.raw / hic_count.hela_g_hic_ActD.total * 300e6
  
  # make inter-dist table
  g_hic_dist_table = make.dist.table(g_hic_matrix.fix,binsize)
  g_hic_ActD_dist_table = make.dist.table(g_hic_ActD_matrix.fix,binsize)
  
  ylim = c(0,7)
  xlim = c(5,8.5)
  
  par(mar=c(2,2,2,1))
  plot(x=log10(g_hic_dist_table$start),y=log10(g_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#1533AD",frame.plot=F,yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  axis(side=2,at=c(0:7))
  
  par(new=T)
  plot(x=log10(g_hic_ActD_dist_table$start),y=log10(g_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#6F81D6",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  title(main = sprintf("HeLa Hi-C chr%d",chrom_index))
  
}


###################################################################
# plot Hela-n-hic & Hela-n-hic-ActD interction-distance figure 
###################################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/make_dist_table.RData")
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")

binsize = 100e3

fig.facet <- layout(matrix(c(1:24),nrow = 6,byrow = FALSE),width = rep(10,24),heights = rep(10,24))
layout.show(fig.facet)

for(chrom_index in c(1:23)){
  print(chrom_index)
  
  load(file=sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_all_matrix_chr_%d_%d.RData",chrom_index,binsize))
  
  # conver as.matrix
  n_hic_matrix.raw = as.matrix(n_hic_matrix)
  n_hic_ActD_matrix.raw = as.matrix(n_hic_ActD_matrix)
  
  # make matrix with the same sequencing depth
  n_hic_matrix.fix = n_hic_matrix.raw / hic_count.hela_n_hic.total * 300e6
  n_hic_ActD_matrix.fix = n_hic_ActD_matrix.raw / hic_count.hela_n_hic_ActD.total * 300e6
  
  # make inter-dist table
  n_hic_dist_table = make.dist.table(n_hic_matrix.fix,binsize)
  n_hic_ActD_dist_table = make.dist.table(n_hic_ActD_matrix.fix,binsize)
  
  ylim = c(0,7)
  xlim = c(5,8.5)
  
  par(mar=c(2,2,2,1))
  plot(x=log10(n_hic_dist_table$start),y=log10(n_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FF8B00",frame.plot=F,yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  par(new=T)
  plot(x=log10(n_hic_ActD_dist_table$start),y=log10(n_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FFC640",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  axis(side=2,at=c(0:7))
  
  title(main = sprintf("HeLa NHi-C chr%d",chrom_index))
}


