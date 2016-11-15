########################################################################################################################
# 绘制Hela-g-hic Hela-n-hic Hela-n-hic/Hela-g-hic matrix 差异图
# 2016-08-31 Howard MENG
########################################################################################################################
rm(list=ls())

###############################################################
# load data & function
###############################################################
load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")
hela_ab_value_table = read.table("~/menghw_HD/data_table/Hela_AB_compart_value.table",header = T,sep = "\t")
hela_all_table = read.table("~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep = "\t")

#################################
# sequencing data
#################################
hela_n_seq_total = 126900190 * 135
hela_n_seq_ActD_total = 282845496 * 135
hela_g_seq_total = 187807841 * 135

hela_g_hic_trans = 55999521 * 2 * 36
hela_g_hic_cis = 247548509 * 2 * 36

hela_n_hic_trans = 62986764 * 2 * 36
hela_n_hic_cis = 273219336 * 2 * 36

hela_n_hic_ActD_trans = 153207862 * 2 * 36
hela_n_hic_ActD_cis = 137755064 * 2 * 36


binsize = 100e3
### 使用coverage进行cutoff
for(chrom_index in c(16:16)){
  chrom_name = chrom.name(chrom_index)
  chrom.region=c(0,chrom.length(chrom_index))
  chrom_start = chrom.region[1]
  chrom_end = chrom.region[2]
  
  image_path = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_%d_%d.RData",chrom_index,binsize)
  load(image_path)
  
  
  filter_vector = hela_all_table$chr_index==chrom_name & hela_all_table$region_min >= chrom.region[1] & hela_all_table$region_max <= chrom.region[2]
  hela_g_seq_count = hela_all_table$g_count[filter_vector]
  hela_g_seq_count.PKM = hela_g_seq_count / (sum(hela_all_table$g_count) / 1e6 * 100) 
  hela_g_seq_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_g_seq_count.PKM)
  
  hela_n_seq_count = hela_all_table$n_count[hela_all_table$chr_index==chrom_name & hela_all_table$region_min >= chrom.region[1] & hela_all_table$region_max <= chrom.region[2]]
  hela_n_seq_count.PKM = hela_n_seq_count / (sum(hela_all_table$n_count) / 1e6 * 100) 
  hela_n_seq_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_n_seq_count.PKM)
  
  hela_g_hic_count = hela_all_table$g_hic[hela_all_table$chr_index==chrom_name & hela_all_table$region_min >= chrom.region[1] & hela_all_table$region_max <= chrom.region[2]]
  hela_g_hic_count.PKM = hela_g_hic_count / (sum(hela_all_table$g_hic) / 1e6 * 100 )
  hela_g_hic_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_g_hic_count.PKM)
  
  hela_n_hic_count = hela_all_table$n_hic[hela_all_table$chr_index==chrom_name & hela_all_table$region_min >= chrom.region[1] & hela_all_table$region_max <= chrom.region[2]]
  hela_n_hic_count.PKM = hela_n_hic_count / (sum(hela_all_table$n_hic) / 1e6 * 100 )
  hela_n_hic_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_n_hic_count.PKM)
  
  hela_n_seq_ratio.bed = hela_n_seq_count.PKM.bed
  hela_n_seq_ratio.bed[,4] = hela_n_seq_count.PKM / hela_g_seq_count.PKM 
  
  hela_n_hic_ratio.bed = hela_n_hic_count.PKM.bed
  hela_n_hic_ratio.bed[,4] = hela_n_hic_count.PKM / hela_g_hic_count.PKM * 1.3
  
  g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,binsize = binsize)
  n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,binsize = binsize)
  
  g_hic_matrix.norm.part.fix = g_hic_matrix.norm.part
  g_hic_matrix.norm.part.fix[g_hic_matrix.norm.part.fix==0] = 1 
  differ_hic_matrix = n_hic_matrix.raw.part / g_hic_matrix.norm.part.fix * 0.95
  differ_hic_matrix[g_hic_matrix.norm.part.fix==0] = 0
  
  png_path = sprintf("~/menghw_HD/R_image/Hela-hic-figure/20160831/Hela-hic_fig1_chr_%d_%d.png",chrom_index,binsize)
  png(png_path,width = 3200,height = 1600)
  
  fig.facet <- layout(matrix(c(1:12),nrow = 4,ncol = 3),width = rep(10,10),heights = rep(c(10,2,2,2),3))
  layout.show(fig.facet)
  
  ###############################
  ## first column
  ###############################
  
  # hela-g-hic matrix
  par(mar=c(0.5,2,0.5,2),family="Arial",font=1)
  plot.matrix(g_hic_matrix.norm,bound.max = 0.95)
  
  # hela-g-hic.norm A/B compartment
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  ab_x = start(pc.norm$PC1)
  ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
  plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
  
  # hela-g-hic coverage PKM
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_g_hic_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#A66500",ylim=c(0,2))
  axis(side=2)
  
  # hela-g-seq coverage PKM
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_g_seq_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#FF9C00",ylim=c(0,2))
  axis(side=2)
  
  
  ###############################
  ## second column
  ###############################
  # hela-g-hic matrix
  par(mar=c(0.5,2,0.5,2),family="Arial",font=1)
  plot.matrix(n_hic_matrix,bound.max = 0.95)
  
  # hela-g-hic.norm A/B compartment
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  ab_x = start(pc.norm$PC1)
  ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
  plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
  
  # hela-g-hic coverage PKM
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_n_hic_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#0E53A7",ylim=c(0,2))
  axis(side=2)
  
  # hela-g-seq coverage PKM
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_n_seq_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#6899D3",ylim=c(0,2))
  axis(side=2)
  
  
  ###############################
  ## third column
  ###############################
  par(mar=c(0.5,2,0.5,2),family="Arial",font=1)
  plot.matrix(differ_hic_matrix,col.min = "navyblue",col.max = "red",bound.max=1,col.boundary = 1.5)
  
  # hela-g-hic.norm A/B compartment
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  ab_x = start(pc.norm$PC1)
  ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
  plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
  
  # ratio
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_n_hic_ratio.bed,chrom_index,chrom_start,chrom_end,track_color = ifelse(hela_n_hic_ratio.bed[,4]>2,"#FF7304","#6899D3"),ylim=c(0,6))
  axis(side=2)
  
  
  # ratio
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_n_seq_ratio.bed,chrom_index,chrom_start,chrom_end,track_color = ifelse(hela_n_hic_ratio.bed[,4]>2.2,"#303030","#6899D3"),ylim=c(0,6))
  axis(side=2)
  
  dev.off()
  
  print(sprintf("The chr%d is done!",chrom_index))
}



########################################################################
### 使用column sum进行cutoff
########################################################################

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_start = 0
  chrom_end = chrom.length(chrom_index)
  chrom.region=c(chrom_start,chrom_end %/% binsize * binsize + binsize )
  
  
  image_path = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_%d_%d.RData",chrom_index,binsize)
  load(image_path)
  
  g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,binsize = binsize)
  n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,binsize = binsize)
  g_hic_matrix.raw.part = matrix.part(g_hic_matrix,chrom_index = chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,binsize = binsize)
  
  g_hic_matrix.norm.part.fix = g_hic_matrix.norm.part
  g_hic_matrix.norm.part.fix[g_hic_matrix.norm.part.fix==0] = 1 
  differ_hic_matrix = n_hic_matrix.raw.part / g_hic_matrix.norm.part.fix * 0.95
  differ_hic_matrix[g_hic_matrix.norm.part.fix==0] = 0
  
  filter_vector = (hela_all_table$chr_index==chrom_name) & (hela_all_table$region_min >= chrom.region[1]) & (hela_all_table$region_max <= chrom.region[2])
  hela_g_seq_count = hela_all_table$g_count[filter_vector]
  hela_g_seq_count.PKM = hela_g_seq_count / (sum(hela_all_table$g_count) / 1e6 * 100) 
  hela_g_seq_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_g_seq_count.PKM)
  
  hela_n_seq_count = hela_all_table$n_count[hela_all_table$chr_index==chrom_name & hela_all_table$region_min >= chrom.region[1] & hela_all_table$region_max <= chrom.region[2]]
  hela_n_seq_count.PKM = hela_n_seq_count / (sum(hela_all_table$n_count) / 1e6 * 100) 
  hela_n_seq_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_n_seq_count.PKM)
  
  hela_g_hic_count = hela_all_table$g_hic[hela_all_table$chr_index==chrom_name & hela_all_table$region_min >= chrom.region[1] & hela_all_table$region_max <= chrom.region[2]]
  hela_g_hic_count.PKM = colSums(g_hic_matrix.raw.part) / (sum(g_hic_matrix.raw.part) / 1e6 * 100 )
  hela_g_hic_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_g_hic_count.PKM)
  
  hela_n_hic_count = hela_all_table$n_hic[hela_all_table$chr_index==chrom_name & hela_all_table$region_min >= chrom.region[1] & hela_all_table$region_max <= chrom.region[2]]
  hela_n_hic_count.PKM = colSums(n_hic_matrix.raw.part) / (sum(n_hic_matrix.raw.part) / 1e6 * 100 )
  hela_n_hic_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_n_hic_count.PKM)
  
  hela_n_seq_ratio.bed = hela_n_seq_count.PKM.bed
  hela_n_seq_ratio.bed[,4] = hela_n_seq_count.PKM / hela_g_seq_count.PKM 
  hela_n_seq_ratio.bed = hela_n_seq_ratio.bed[which(is.finite(hela_n_seq_ratio.bed[,4])),]
  
  hela_n_hic_ratio.bed = hela_n_hic_count.PKM.bed
  hela_n_hic_ratio.bed[,4] = hela_n_hic_count.PKM / hela_g_hic_count.PKM
  hela_n_hic_ratio.bed = hela_n_hic_ratio.bed[which(is.finite(hela_n_hic_ratio.bed[,4])),]
  

  png_path = sprintf("~/menghw_HD/R_image/Hela-hic-figure/20160831/Hela-hic_fig1_chr_%d_%d_v2.png",chrom_index,binsize)
  png(png_path,width = 3200,height = 1600)
  
  fig.facet <- layout(matrix(c(1:12),nrow = 4,ncol = 3),width = rep(10,10),heights = rep(c(10,2,2,2),3))
  layout.show(fig.facet)
  
  ###############################
  ## first column
  ###############################
  
  # hela-g-hic matrix
  par(mar=c(0.5,2,0.5,2),family="Arial",font=1)
  plot.matrix(g_hic_matrix.norm,bound.max = 0.95)
  
  # hela-g-hic.norm A/B compartment
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  ab_x = start(pc.norm$PC1)
  ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
  plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
  
  # hela-g-hic coverage PKM
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_g_hic_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#A66500",ylim=c(0,40))
  axis(side=2)
  
  # hela-g-seq coverage PKM
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_g_seq_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#FF9C00",ylim=c(0,2))
  axis(side=2)
  
  
  ###############################
  ## second column
  ###############################
  # hela-g-hic matrix
  par(mar=c(0.5,2,0.5,2),family="Arial",font=1)
  plot.matrix(n_hic_matrix,bound.max = 0.95)
  
  # hela-g-hic.norm A/B compartment
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  ab_x = start(pc.norm$PC1)
  ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
  plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
  
  # hela-g-hic coverage PKM
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_n_hic_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#0E53A7",ylim=c(0,40))
  axis(side=2)
  
  # hela-g-seq coverage PKM
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_n_seq_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#6899D3",ylim=c(0,2))
  axis(side=2)
  
  
  ###############################
  ## third column
  ###############################
  NAT_cutoff = 1.5
  NAD_cutoff = 2
  
  par(mar=c(0.5,2,0.5,2),family="Arial",font=1)
  plot.matrix(differ_hic_matrix,col.min = "navyblue",col.max = "red",bound.max=1,col.boundary = 2)
  
  # hela-g-hic.norm A/B compartment
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  ab_x = start(pc.norm$PC1)
  ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
  plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
  
  # ratio
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_n_hic_ratio.bed,chrom_index,chrom_start,chrom_end,track_color = ifelse(hela_n_hic_ratio.bed[,4]>NAT_cutoff,"#FF7304","#6899D3"),ylim=c(0,6))
  axis(side=2)
  
  
  # ratio
  par(mar=c(0.5,2,0,2),family="Arial",font=1)
  plot.barplot(hela_n_seq_ratio.bed,chrom_index,chrom_start,chrom_end,track_color = ifelse(hela_n_hic_ratio.bed[,4]>NAD_cutoff,"#303030","#6899D3"),ylim=c(0,6))
  axis(side=2)
  
  dev.off()
  
  print(sprintf("The chr%d is done!",chrom_index))
}





###############################################################
# plot 4 types matrix 
###############################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")

binsize = 100e3

for(chrom_index in c(1:23)){
  
  png_path = sprintf("~/menghw_HD/R_image/Hela-hic-figure/20160908/Hela-hic_matrix_chr_%d_%d_v1.png",chrom_index,binsize)
  png(png_path,width = 2000,height = 2000)
  
  chrom_name = chrom.name(chrom_index)
  image_path = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_%d_%d.RData",chrom_index,binsize)
  load(image_path)
  
  matrix_path = sprintf("~/menghw_HD/our_data/Hela-genome-hic-actD-hg19/matrix/raw_matrix/chr_%d_%d_MAPQ20.txt",chrom_index,binsize)
  g_hic_ActD_matrix = as.matrix(read.table(matrix_path,header = F,sep = ","))
  
  fig.facet <- layout(matrix(c(1:4),nrow = 2,ncol = 2),width = rep(10,4),heights = rep(10,4))
  layout.show(fig.facet)
  
  par(mar=c(1,1,1,1),family="Arial",font=1)
  plot.matrix(g_hic_matrix.norm,bound.max = 0.95)
  
  par(mar=c(1,1,1,1),family="Arial",font=1)
  plot.matrix(g_hic_ActD_matrix,bound.max = 0.95)
  
  par(mar=c(1,1,1,1),family="Arial",font=1)
  plot.matrix(n_hic_matrix,bound.max = 0.95)
  
  par(mar=c(1,1,1,1),family="Arial",font=1)
  plot.matrix(n_hic_ActD_matrix,bound.max = 0.95)
  
  dev.off()
  print(sprintf("The chr%d is done!",chrom_index))
}







