########################################################################################################################
# plot B-compartment among chromosomes
# 2016-10-31 Howard MENG
########################################################################################################################
rm(list=ls())

###############################################################
# load data & function
###############################################################
load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")
hela_ab_value_table = read.table("~/menghw_HD/data_table/Hela_AB_compart_value.table",header = T,sep = "\t")
hela_all_table = read.table("~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep = "\t")
hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header=F,sep="\t")

binsize = 100e3

png(file="~/menghw_HD/R_image/Hela-hic-figure/Hela_NAIR_B_compartment.png",width = 1600,height = 2400)

fig.facet <- layout(matrix(c(1:96),nrow = 96,ncol = 1),width = rep(10,96),heights = rep(c(0.3,1,0.5,0.5),24))
layout.show(fig.facet)


for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  print(chrom_name)
  chrom.region=c(0,chrom.length(chrom_index))
  chrom_start = chrom.region[1]
  chrom_end = chrom.region[2]
  
  image_path = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_%d_%d.RData",chrom_index,binsize)
  load(image_path)
  
  if(chrom_index > 12){
    par.mar.left = 3
  }else{
    par.mar.left = 3
  }
  
  # hg19 cytoband
  par(mar=c(0,par.mar.left,0.5,0.5))
  plot.chromosome(chrom_index,0,chrom.length(1))

  # hela-g-hic.norm A/B compartment
  par(mar=c(0,par.mar.left,0.5,0.5),family="Arial",font=1)
  ab_x = start(pc.norm$PC1)
  ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
  plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(0,chrom.length(1)),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
  
  # hela NAIR
  par(mar=c(1,par.mar.left,0,0.5),family="Arial",font=1)
  plot.peak(df.peak = hela_NAIR_table,chrom_index,chrom_len,chrom_start = 0,chrom_end =chrom.length(1) ,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  par(mar=c(0,par.mar.left,0,0),family="Arial",font=1)
  plot(x=1,y=1,frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
  
}

dev.off()






