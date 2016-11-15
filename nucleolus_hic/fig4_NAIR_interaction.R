######################################################
# NAIR 区的region plot inter-dist  
######################################################
rm(list=ls())
hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
colnames(hela_NAIR_table) = c("chrom_name","start","end","peak_index","info")

# fix function 
make.NAIR.dist.table <- function(hic_matrix,binsize,NAIR_table,chrom_name){
  # 输入chromosome hic matrix 与binsize 生成distance-interaction dataframe
  chrom_NAIR_table = NAIR_table[NAIR_table[,1]==chrom_name,]
  colnames(chrom_NAIR_table) = c("chrom_name","start","end","peak_index","info")
  NAIR_bin_index_list = NULL
  for(i in c(1:nrow(chrom_NAIR_table))){
    bin.index.start = chrom_NAIR_table$start[i] %/% binsize 
    bin.index.end = chrom_NAIR_table$end[i] %/% binsize
    NAIR_bin_index_list = c(NAIR_bin_index_list,c(bin.index.start:bin.index.end))
  }
  
  mat.nrow = dim(hic_matrix)[1]
  mat.ncol = dim(hic_matrix)[2]
  
  count_vector = rep(0,mat.nrow)
  for(row.index in c(1:mat.nrow)){
    x = c(row.index : mat.nrow)
    y = c(1:length(x))
    inter_count = 0
    for(xy.index in c(1:length(x))){
      x0 = x[xy.index]
      y0 = y[xy.index]
      if(x0 %in% NAIR_bin_index_list & y0 %in% NAIR_bin_index_list){
        inter_count = inter_count + hic_matrix[x0,y0]  
      }
    }
    count_vector[row.index] = inter_count
  }
  
  # 生成返回的data frame
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }else if(chrom_index == 25){
    chrom_name = "rDNA"
  }
  
  start_vector = as.integer(c(1:mat.nrow) * binsize)
  
  inter_dist_table = data.frame(chrom_name = rep(chrom_name,length(count_vector)),
                                start = start_vector,
                                hic_count = count_vector)
  return(inter_dist_table)
}


###### plot 
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")
load(file="~/menghw_HD/R_code/my_function/make_dist_table.RData")
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")
hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
colnames(hela_NAIR_table) = c("chrom_name","start","end","peak_index","info")
binsize = 100e3

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  print(chrom_name)
  
  load(file=sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_all_matrix_chr_%d_%d.RData",chrom_index,binsize))
  
  # conver as.matrix
  g_hic_matrix.raw = as.matrix(g_hic_matrix)
  n_hic_matrix.raw = as.matrix(n_hic_matrix)
  g_hic_ActD_matrix.raw = as.matrix(g_hic_ActD_matrix)
  n_hic_ActD_matrix.raw = as.matrix(n_hic_ActD_matrix)
  
  # make matrix with the same sequencing depth
  g_hic_matrix.fix = g_hic_matrix.raw / hic_count.hela_g_hic.total * 300e6
  n_hic_matrix.fix = n_hic_matrix.raw / hic_count.hela_n_hic.total * 300e6
  g_hic_ActD_matrix.fix = g_hic_ActD_matrix.raw / hic_count.hela_g_hic_ActD.total * 300e6
  n_hic_ActD_matrix.fix = n_hic_ActD_matrix.raw / hic_count.hela_n_hic_ActD.total * 300e6
  
  # make inter-dist table
  print("g_hic")
  g_hic_dist_table = make.NAIR.dist.table(g_hic_matrix.fix,binsize,hela_NAIR_table,chrom_name)
  print("g_hic_ActD")
  g_hic_ActD_dist_table = make.NAIR.dist.table(g_hic_ActD_matrix.fix,binsize,hela_NAIR_table,chrom_name)
  print("n_hic")
  n_hic_dist_table = make.NAIR.dist.table(n_hic_matrix.fix,binsize,hela_NAIR_table,chrom_name)
  print("n_hic_ActD")
  n_hic_ActD_dist_table = make.NAIR.dist.table(n_hic_ActD_matrix.fix,binsize,hela_NAIR_table,chrom_name)
  
  # save image
  save.image(file = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela_inter_dist/hela_inter_dist_chr_%d_%d.RData",chrom_index,binsize))
  print(sprintf("chr%d work is done!",chrom_index))
  
#   ylim = c(0,7)
#   xlim = c(0,chrom.length(chrom_index))
#   
#   png(file=sprintf("~/menghw_HD/R_image/Hela-hic-figure/inter_dist_plot/chr_%d_%d_NAIR_dist_v2.png",chrom_index,binsize),width = 1000,height = 1000)
#   
#   plot(x=g_hic_dist_table$start,y=log10(g_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#1533AD",frame.plot=F,yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
#   axis(side=2,at=c(0:7))
#   
#   par(new=T)
#   plot(x=g_hic_ActD_dist_table$start,y=log10(g_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#6F81D6",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
#   
#   par(new=T)
#   plot(x=n_hic_dist_table$start,y=log10(n_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FF8B00",frame.plot=F,yaxt="n",xlab="",xaxt="n",ylab="",xlim=xlim,ylim=ylim)
#   
#   par(new=T)
#   plot(x=n_hic_ActD_dist_table$start,y=log10(n_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FFC640",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
#   
#   dev.off()
  
}
