##################################################################
# use png test B compartment and enrichment region
##################################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")
# hela_NAD_table = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep = "\t")
# hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
hela_ab_value_table = read.table("~/menghw_HD/data_table/Hela_AB_compart_value.table",header = T,sep = "\t")

binsize = 100e3
for(chrom_index in c(10:21)){
  chrom_name = chrom.name(chrom_index)
  chrom_start = 0 
  chrom_end = chrom.length(chrom_index)
  
  image_path = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_track_chr_%d_%d.RData",chrom_index,binsize)
  load(image_path)
  png_path = sprintf("~/menghw_HD/R_image/Hela-hic-figure/20160830/Hela_matrix_chr_%d_%d_v1.png",chrom_index,binsize)
  
  png(png_path,width = 1400,height = 3000)
  fig.facet <- layout(matrix(c(1:6),nrow = 6,byrow = FALSE),width = rep(10,6),heights = c(10,1,1,1,10,1))
  layout.show(fig.facet)
  
  ## g_hic_matrix.ice.norm
  par(mar=c(0.5,1,1,1))
  plot.matrix(matrix.part(g_hic_matrix.norm,chrom_index,chrom_start,chrom_end,binsize),max_bound = 0.96)
  # plot.matrix(g_hic_matrix.norm,max_bound = 0.95)
  # box(lwd=3)
  
  ## A,B compartment
  par(mar=c(0.5,1,0,1))
  ab_x = start(pc.norm$PC1)
  ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
  plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",yaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
  
  ## NAD
  par(mar=c(0.5,1,0,1))
  plot.peak(hela_NAD_table,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#EE5C42",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
  
  ## NAT
  par(mar=c(0.5,1,0,1))
  plot.peak(hela_NAT_table,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#EE5C42",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
  
  ##n_hic_matrix
  par(mar=c(0.1,1,0,1))
  plot.matrix(matrix.part(n_hic_matrix,chrom_index,chrom_start,chrom_end,binsize),max_bound = 0.96)
  
  ## chrom 坐标轴
  par(mar=c(0.5,1,0.5,1))
  plot.chromosome(chrom_index,chrom_start,chrom_end,border.lwd = 2)
  
  dev.off()
  print(sprintf("The chr%d ice.norm figure is done!",chrom_index))
}

