########################################################################################################################
# 绘制Hela-g-hic Hela-n-hic Hela-n-hic/Hela-g-hic matrix 差异图
# 2016-09-30 Howard MENG
########################################################################################################################

###########################################################################################################
# fig1 B identify NAIRs 
###########################################################################################################

rm(list=ls())

###############################################################
# load data & function
###############################################################
load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")
hela_ab_value_table = read.table("~/menghw_HD/data_table/Hela_AB_compart_value.table",header = T,sep = "\t")
hela_all_table = read.table("~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep = "\t")


binsize = 100e3
### 使用coverage进行cutoff
for(chrom_index in c(16:16)){
  chrom_name = chrom.name(chrom_index)
  # chrom.region=c(0,chrom.length(chrom_index))
  chrom.region=c(50e6,90e6)
  chrom_start = chrom.region[1]
  chrom_end = chrom.region[2]

  load(sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_%d_%d.RData",chrom_index,binsize))
  
  
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
  g_hic_matrix.raw.part = matrix.part(g_hic_matrix,chrom_index = chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,binsize = binsize)
  n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,binsize = binsize)
  
  g_hic_matrix.norm.part.fix = g_hic_matrix.norm.part
  g_hic_matrix.norm.part.fix[g_hic_matrix.norm.part.fix==0] = 1 
  differ_hic_matrix = n_hic_matrix.raw.part / g_hic_matrix.norm.part.fix * 0.95
  differ_hic_matrix[g_hic_matrix.norm.part.fix==0] = 0
  differ_hic_matrix.part = n_hic_matrix.raw.part / g_hic_matrix.raw.part
  differ_hic_matrix.part[n_hic_matrix.raw.part==0 | g_hic_matrix.raw.part==0] = 0
  
  png_path = sprintf("~/menghw_HD/R_image/Hela-hic-figure/Hela-hic_fig1B.png")
  png(png_path,width = 1600,height  = 800)
  
  fig.facet <- layout(matrix(c(1:12),nrow = 3,ncol = 3),width = rep(10,10),heights = rep(c(10,2,2),3))
  layout.show(fig.facet)
  
  ###############################
  ## first column
  ###############################
  
  # hela-g-hic matrix
  par(mar=c(0.5,2,0.5,2),family="Arial",font=1)
  plot.matrix(g_hic_matrix.norm.part,bound.max = 0.95)
  
  # # hela-g-hic.norm A/B compartment
  # par(mar=c(0.5,2,0,2),family="Arial",font=1)
  # ab_x = start(pc.norm$PC1)
  # ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
  # plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
  
  # hela-g-hic coverage PKM
  par(mar=c(0.5,3,0,0),family="Arial",font=1)
  plot.barplot(hela_g_hic_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#A66500",ylim=c(0,2))
  axis(side=2,at=c(0,1,2),labels = F,lwd = 3,tck=-0.1)
  
  # hela-g-seq coverage PKM
  par(mar=c(0.5,3,0,0),family="Arial",font=1)
  plot.barplot(hela_g_seq_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#FF9C00",ylim=c(0,2))
  axis(side=2,at=c(0,1,2),labels = F,lwd = 3,tck=-0.1)
  # axis(side=2)
  
  
  ###############################
  ## second column
  ###############################
  # hela-g-hic matrix
  # par(mar=c(0.5,2,0.5,2),family="Arial",font=1)
  par(mar=c(1,1,1,1),family="Arial",font=1)
  plot.matrix(n_hic_matrix.raw.part,bound.max = 0.95)
  
  # hela-g-hic.norm A/B compartment
  # par(mar=c(0.5,2,0,2),family="Arial",font=1)
  # ab_x = start(pc.norm$PC1)
  # ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
  # plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
  
  # hela-g-hic coverage PKM
  # par(mar=c(0.5,2,0,2),family="Arial",font=1)
  par(mar=c(0.5,3,0,0),family="Arial",font=1)
  plot.barplot(hela_n_hic_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#0E53A7",ylim=c(0,2))
  axis(side=2,at=c(0,1,2),labels = F,lwd = 3,tck=-0.1)
  # axis(side=2)
  
  # hela-g-seq coverage PKM
  # par(mar=c(0.5,2,0,2),family="Arial",font=1)
  par(mar=c(0.5,3,0,0),family="Arial",font=1)
  plot.barplot(hela_n_seq_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#6899D3",ylim=c(0,2))
  axis(side=2,at=c(0,1,2),labels = F,lwd = 3,tck=-0.1)
  # axis(side=2)
  
  
  ###############################
  ## third column
  ###############################
  # par(mar=c(0.5,2,0.5,2),family="Arial",font=1)
  par(mar=c(1,1,1,1),family="Arial",font=1)
  plot.matrix(differ_hic_matrix.part,col.min = "navyblue",col.max = "red",bound.max=0.99,col.boundary = 1.414)
  
  differ_hic_matrix.part.log = log2(differ_hic_matrix.part)
  differ_hic_matrix.part.log[differ_hic_matrix.part==0] = -3
  quantile(differ_hic_matrix.part.log,prob=seq(0,1,0.01))
  plot.matrix(differ_hic_matrix.part.log,col.min = "navyblue",col.max = "red",bound.min = 0.25,bound.max=0.99,col.boundary = 0)
  
  # hela-g-hic.norm A/B compartment
  # par(mar=c(0.5,2,0,2),family="Arial",font=1)
  # ab_x = start(pc.norm$PC1)
  # ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
  # plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
  
  # ratio
  # par(mar=c(0.5,2,0,2),family="Arial",font=1)
  par(mar=c(0.5,3,0,0),family="Arial",font=1)
  plot.barplot(hela_n_hic_ratio.bed,chrom_index,chrom_start,chrom_end,track_color = ifelse(hela_n_hic_ratio.bed[,4]>2,"#FF7304","#6899D3"),ylim=c(0,4))
  axis(side=2,at=c(0,2,4),labels = F,lwd = 3,tck=-0.1)
  
  # ratio
  # par(mar=c(0.5,2,0,2),family="Arial",font=1)
  par(mar=c(0.5,3,0,0),family="Arial",font=1)
  plot.barplot(hela_n_seq_ratio.bed,chrom_index,chrom_start,chrom_end,track_color = ifelse(hela_n_hic_ratio.bed[,4]>2.2,"#303030","#6899D3"),ylim=c(0,4))
  axis(side=2,at=c(0,2,4),labels = F,lwd = 3,tck=-0.1)
  # axis(side=2)
  
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



###########################################################################################################
# fig1 E data correlation
###########################################################################################################
rm(list=ls())
hela_all_table = read.table("~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep = "\t")
relation_matrix = matrix(rep(1,16),4)
relation_matrix[1,2] = cor(x=hela_all_table$g_count,y=hela_all_table$g_hic)
relation_matrix[2,1] = relation_matrix[1,2]

relation_matrix[1,3] = cor(x=hela_all_table$g_count,y=hela_all_table$n_count)
relation_matrix[3,1] = relation_matrix[1,3]

relation_matrix[1,4] = cor(x=hela_all_table$g_count,y=hela_all_table$n_hic)
relation_matrix[4,1] = relation_matrix[1,4]

relation_matrix[2,3] = cor(x=hela_all_table$g_hic,y=hela_all_table$n_count)
relation_matrix[3,2] = relation_matrix[2,3]


relation_matrix[2,4] = cor(x=hela_all_table$g_hic,y=hela_all_table$n_hic)
relation_matrix[4,2] = relation_matrix[2,4]

relation_matrix[3,4] = cor(x=hela_all_table$n_count,y=hela_all_table$n_hic)
relation_matrix[4,3] = relation_matrix[3,4]

# colnames(relation_matrix) = c("WGS","Hi-C","NS","nHi-C")
# rownames(relation_matrix) = colnames(relation_matrix)
library(pheatmap)
library(gplots)
pheatmap(relation_matrix,color = colorpanel(100,low = "yellow",high = "red"),cluster_rows = F,cluster_cols = F,cex=1.5,display_numbers=T,fontsize_number = 18,number_color = "black")


###########################################################################################################
# fig1 F NAIR interaction 
###########################################################################################################
barplot.mat = c(0.3447659,0.6552341,0.6072660,0.3927340)
barplot.mat = matrix(barplot.mat,2,2)
# barplot.mat[1,1] = g_hic_info_NAT[2] / g_hic_info_NAT[1]
# barplot.mat[2,1] = g_hic_info_NAT[3] / g_hic_info_NAT[1]

# barplot.mat[1,2] = n_hic_info_NAT[2] / n_hic_info_NAT[1]
# barplot.mat[2,2] = n_hic_info_NAT[3] / n_hic_info_NAT[1]

####################################
## NAD-NAD interaction barplot 
####################################
par(mar=c(2,4,1,1),family="Arial",font=1)
barplot(barplot.mat,axes = F,col=c("#FF7304","#7080D7"),width = 0.5,space = 0.25,border = F)
axis(side=2,at=seq(0,1,0.25),labels = F,cex.axis=3,lwd=3,tck=-0.05)
# text(x=-0.1,y=seq(0,1,0.25),labels = paste0(seq(0,100,25),"%"),cex = 3,xpd=T,pos = NULL)
# text(x=c(0.45,1.05),y=-0.05,labels = c("In Situ Hi-C","N-Hi-C"),cex = 3,xpd=T,srt = 45,pos = 2)







