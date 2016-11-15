##################################################################################
# our NAD and older NAD overlap
##################################################################################
rm(list=ls())
our_NAD_table = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep="\t")
old_NAD_table = read.table("~/menghw_HD/data_table/2010_NAD_region/2010_NAD_hg19_w1e5.bed",header = F,sep = "\t")

sum(our_NAD_table$V3 - our_NAD_table$V2)
sum(old_NAD_table$V3 - old_NAD_table$V2)

load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")
chrom_name_list = c(paste("chr",c(1:22),sep = ""),"chrX")
NAD_overlap_table = overlap(our_NAD_table,old_NAD_table,level_vetor = chrom_name_list)

sum(NAD_overlap_table$end - NAD_overlap_table$start)


png(filename = "~/menghw_HD/R_image/Hela-hic-figure/NAD_NAD_overlap.png",width = 2000,height = 3000)

fig.facet <- layout(matrix(c(1:92),nrow = 92,byrow = FALSE),width = rep(10,92),heights = rep(c(1,1.5,1,1),23))
layout.show(fig.facet)

chrom_start = 0
chrom_end = chrom.length(1)

for(chrom_index in c(1:23)){
  ## our NAD
  par(mar=c(0.5,1,0,1))
  plot.peak(our_NAD_table,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#472B83",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
  
  ## chrom 坐标轴
  par(mar=c(0.5,1,0,1))
  plot.chromosome(chrom_index,chrom_start,chrom_end,border.lwd = 2)
  
  ## old NAD
  par(mar=c(0.5,1,0,1))
  plot.peak(old_NAD_table,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#006957",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
  
  
  par(mar=c(0.5,1,0,1))
  plot(x=c(chrom_start,chrom_end),y=c(-1,1),type="n",frame.plot = F,ylab = "",xlab = "",xaxt = "n",yaxt="n")

}


dev.off()



##################################################################################
# Hela-NAIR v.s. U2OS-NAIR
##################################################################################
rm(list=ls())
hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep="\t")
u2os_NAIR_table = read.table("~/menghw_HD/data_table/u2os_NAIR_data/U2OS_NAIR_v1_fix.bed",header = T,sep = "\t")
colnames(hela_NAIR_table) = colnames(u2os_NAIR_table)


load(file="~/menghw_HD/R_code/my_function/MyPlot_v07.RData")
chrom_name_list = c(paste("chr",c(1:22),sep = ""),"chrX")
NAIR_overlap_table = overlap(hela_NAIR_table,u2os_NAIR_table,level_vetor = chrom_name_list)

sum(NAIR_overlap_table$end - NAIR_overlap_table$start) / sum(hela_NAIR_table$end - hela_NAIR_table$start)
sum(NAIR_overlap_table$end - NAIR_overlap_table$start) / sum(u2os_NAIR_table$end - u2os_NAIR_table$start)

png(filename = "~/menghw_HD/R_image/Hela-hic-figure/NAIR_NAIR_overlap.png",width = 2000,height = 3000)

fig.facet <- layout(matrix(c(1:92),nrow = 92,byrow = FALSE),width = rep(10,92),heights = rep(c(1,1.5,1,1),23))
layout.show(fig.facet)

chrom_start = 0
chrom_end = chrom.length(1)

for(chrom_index in c(1:23)){
  ## Hela NAIR
  par(mar=c(0.5,1,0,1))
  plot.peak(hela_NAIR_table,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
  
  ## chrom 坐标轴
  par(mar=c(0.5,1,0,1))
  plot.chromosome(chrom_index,chrom_start,chrom_end,border.lwd = 2)
  
  ## U2OS NAIR
  par(mar=c(0.5,1,0,1))
  plot.peak(u2os_NAIR_table,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#A65A00",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
  
  
  par(mar=c(0.5,1,0,1))
  plot(x=c(chrom_start,chrom_end),y=c(-1,1),type="n",frame.plot = F,ylab = "",xlab = "",xaxt = "n",yaxt="n")
  
}


dev.off()


##################################################################################
# NAIR and cLAD overlap
##################################################################################
cLAD_table = read.table("~/menghw_HD/data_table/LAD_region/GSE22428_LAD_overlap_w1e5.txt",header = F,sep = "\t")
colnames(cLAD_table) = c("chrom_name","start","end","peak_index","info")

hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep="\t")
colnames(hela_NAIR_table) = c("chrom_name","start","end","peak_index","info")

load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")

chrom_name_list = c(paste("chr",c(1:22),sep = ""),"chrX")
NAIR_overlap_table = overlap(hela_NAIR_table,cLAD_table,level_vetor = chrom_name_list)


sum(NAIR_overlap_table$end - NAIR_overlap_table$start) / sum(hela_NAIR_table$end - hela_NAIR_table$start)

png(filename = "~/menghw_HD/R_image/Hela-hic-figure/NAIR_cLAD_overlap.png",width = 2000,height = 3000)

fig.facet <- layout(matrix(c(1:92),nrow = 92,byrow = FALSE),width = rep(10,92),heights = rep(c(1,1.5,1,1),23))
layout.show(fig.facet)

chrom_start = 0
chrom_end = chrom.length(1)

for(chrom_index in c(1:23)){
  ## Hela NAIR
  par(mar=c(0.5,1,0,1))
  plot.peak(hela_NAIR_table,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
  
  ## chrom 坐标轴
  par(mar=c(0.5,1,0,1))
  plot.chromosome(chrom_index,chrom_start,chrom_end,border.lwd = 2)
  
  ## cLAD
  par(mar=c(0.5,1,0,1))
  plot.peak(cLAD_table,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#805E15",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
  
  
  par(mar=c(0.5,1,0,1))
  plot(x=c(chrom_start,chrom_end),y=c(-1,1),type="n",frame.plot = F,ylab = "",xlab = "",xaxt = "n",yaxt="n")
  
}


dev.off()
















