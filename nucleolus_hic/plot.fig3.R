#############################################################################################################
# 2016-09-10 Howard MENG
# plot fig3 figures compare NAIR and B compartment
#############################################################################################################

########################################################################
# fig3-A overlap between NAIR and B-compartment
########################################################################
rm(list=ls())

# load data sequencing depth

load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")
hela_ab_value_table = read.table("~/menghw_HD/data_table/Hela_AB_compart_value.table",header = T,sep = "\t")
hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")


# load matrix
chrom_index = 18
# chrom_index = 19
binsize = 50e3
# binsize = 100e3


# image_path = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_all_matrix_chr_%d_%d.RData",chrom_index,binsize)
image_path = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_%d_%d.RData",chrom_index,binsize)
load(image_path)
load(file="~/menghw_HD/R_code/my_function/MyPlot_v07.RData")
chrom_name = chrom.name.factor(chrom_index)

###########################
# get whole matrix first
###########################
segment.col = c(0,chrom.length(chrom_index))
segment.row = segment.col
chrom_start = segment.col[1]
chrom_end = segment.col[2]

###########################
# get part matrix 
###########################
segment.row = c(25e6,40e6)
segment.col = c(45e6,75e6)

segment.row = c(52.5e6,65e6)
segment.col = c(52.5e6,65e6)

segment.row = c(30e6,59e6)
segment.col = c(30e6,59e6)

chrom_start = segment.col[1]
chrom_end = segment.col[2]


###########################
# fix matrix
###########################
g_hic_matrix.part = matrix.part.corner(g_hic_matrix,chrom_index,segment.row,segment.col,binsize)
n_hic_matrix.part = matrix.part.corner(n_hic_matrix,chrom_index,segment.row,segment.col,binsize)

g_hic_matrix.ActD = read.table("~/menghw_HD/our_data/Hela-genome-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_100000/chr_19_100000_MAPQ20.txt",header = F,sep = ",")

g_hic_matrix.ActD.part = matrix.part.corner(g_hic_matrix.ActD,chrom_index,segment.row,segment.col,binsize)
n_hic_matrix.ActD.part = matrix.part.corner(n_hic_ActD_matrix,chrom_index,segment.row,segment.col,binsize)

# g_hic_matrix.part.fix = g_hic_matrix.part
# n_hic_matrix.part.fix = n_hic_matrix.part

g_hic_matrix.part.fix = g_hic_matrix.part / hic_count.hela_g_hic.total * 300e6
n_hic_matrix.part.fix = n_hic_matrix.part / hic_count.hela_n_hic.total * 300e6
g_hic_matrix.ActD.part.fix = g_hic_matrix.ActD.part / hic_count.hela_g_hic_ActD.total * 300e6
n_hic_matrix.ActD.part.fix = n_hic_matrix.ActD.part / hic_count.hela_n_hic_ActD.total * 300e6

diff_matrix.part = n_hic_matrix.part.fix / g_hic_matrix.part.fix
diff_matrix.part[g_hic_matrix.part.fix==0]=0


g_hic_diff_matrix.part = g_hic_matrix.ActD.part.fix / g_hic_matrix.part.fix
g_hic_diff_matrix.part[g_hic_matrix.part.fix==0] = 0

n_hic_diff_matrix.part = n_hic_matrix.ActD.part.fix / n_hic_matrix.part.fix
n_hic_diff_matrix.part[n_hic_matrix.part.fix==0] = 0

###########################
# plot heatmap 
###########################
segment.row.part = c(25e6,40e6)
segment.col.part = c(45e6,75e6)
segment.col.fix = segment.col.part %/% binsize
segment.row.fix = segment.row.part %/% binsize

# png(file="~/menghw_HD/R_image/Hela-hic-figure/20160910/fig4A.png",width = 1200,height = 1200)
fig.facet <- layout(matrix(c(1:4),nrow = 4,byrow = FALSE),width = rep(20,4),heights = c(8,1,1,8))
layout.show(fig.facet)

# hela-g-hic heatmap
par(mar=c(0.5,1,0.5,1),family="Arial",font=1)
# par(mar=c(1,1,1,1),family="Arial",font=1)
# plot.matrix(g_hic_matrix.part.fix,bound.max = 0.95)
plot.matrix(diff_matrix.part,bound.max = 0.90,col.boundary = 1,col.min = "blue")
# plot.matrix(log2(g_hic_matrix.part.fix+1),bound.max = 0.95)
# plot.matrix(log2(n_hic_matrix.part.fix+1),bound.max = 0.995)
segment.row.part.plot = c(25e6,40e6)
segment.col.part.plot = c(45e6,75e6)
segment.col.fix = segment.col.part.plot %/% binsize
segment.row.fix = segment.row.part.plot %/% binsize
plot.segment(segment.row.fix,segment.col.fix,color = "#346A7F")

segment.row.part.plot = c(50e6,65e6)
segment.col.part.plot = c(50e6,65e6)
segment.col.fix = segment.col.part.plot %/% binsize
segment.row.fix = segment.row.part.plot %/% binsize
plot.segment(segment.row.fix,segment.col.fix,color = "#346A7F")


# A/B compartment
par(mar=c(0.5,1,0.5,1),family="Arial",font=1)
ab_x = start(pc.norm$PC1)
ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])

ab_x.fix = ab_x[ab_x >= chrom_start & ab_x <= chrom_end]
ab_y.fix = ab_y[ab_x >= chrom_start & ab_x <= chrom_end]
plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=3,col=ifelse(ab_y>0,"orange","blue"))

# NAIR 
par(mar=c(0.5,1,0.5,1),family="Arial",font=1)
plot.peak(df.peak = hela_NAT_table,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

# hela-n-hic heatmap
plot.matrix(n_hic_matrix.part.fix,bound.max = 0.92)

segment.row.part.plot = c(25e6,40e6)
segment.col.part.plot = c(45e6,75e6)
segment.col.fix = segment.col.part.plot %/% binsize
segment.row.fix = segment.row.part.plot %/% binsize
plot.segment(segment.row.fix,segment.col.fix,color = "#346A7F")

segment.row.part.plot = c(50e6,65e6)
segment.col.part.plot = c(50e6,65e6)
segment.col.fix = segment.col.part.plot %/% binsize
segment.row.fix = segment.row.part.plot %/% binsize
plot.segment(segment.row.fix,segment.col.fix,color = "#346A7F")

dev.off()

# chromosome 18 cytoband
par(mar=c(1,1,1,1))
plot.chromosome(chrom_index,chrom_start,chrom_end,border.lwd = 3)


######################################################
# plot heatmap for g_hic_differ and n_hic_differ
######################################################
fig.facet <- layout(matrix(c(1:2),nrow = 2,byrow = FALSE),width = rep(20,2),heights = c(8,1))
layout.show(fig.facet)

# hela-g-hic heatmap
par(mar=c(0.5,1,0.5,1),family="Arial",font=1)

g_hic_diff_matrix.part.log = log2(g_hic_diff_matrix.part)
g_hic_diff_matrix.part.log[g_hic_diff_matrix.part==0] = 0

n_hic_diff_matrix.part.log = log2(n_hic_diff_matrix.part)
n_hic_diff_matrix.part.log[n_hic_diff_matrix.part==0] = 0

# plot.matrix(g_hic_diff_matrix.part,bound.max = 0.90,col.boundary = 1,col.min = "blue")
plot.matrix(g_hic_diff_matrix.part.log,bound.max = 0.98,bound.min = 0.05,col.boundary = 0,col.min = "blue")
quantile(g_hic_diff_matrix.part.log,0.05)
quantile(g_hic_diff_matrix.part.log,0.98)

plot.matrix(n_hic_diff_matrix.part.log,bound.max = 0.99,bound.min = 0.05,col.boundary = 0,col.min = "blue")
quantile(n_hic_diff_matrix.part.log,0.01)
quantile(n_hic_diff_matrix.part.log,0.99)

# NAIR 
par(mar=c(0.5,1,0.5,1),family="Arial",font=1)
hela_NAT_table.fix = hela_NAT_table[hela_NAT_table$V3 - hela_NAT_table$V2 > 500e3,]
plot.peak(df.peak = hela_NAT_table.fix,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
dev.off()


########################################################################
# fig3 function test plot segment line
########################################################################
plot.matrix(matrix(c(1:100),10))
plot.segment(segment.row = c(5,8),segment.col = c(6,10))


#############################################################################################
# fig3 chromsome NAIR B-compartment 
#############################################################################################
rm(list=ls())

##################################
# loading data 
##################################
# sequencing depth
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")

# A/B compartment with 100Kb binsize
hela_AB_table = read.table("~/menghw_HD/data_table/Hela-genome-hic_AB_all_100000_ice.table",header = T,sep = "\t")

# NAD/NAT region 
hela_NAD_table = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep = "\t")
hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")

# Hi-C column Sum 20Kb binsize
hic_colsum_table = read.table("~/menghw_HD/data_table/Hela_hic_colsum.table",header = T,sep = "\t")

# function 
load(file="~/menghw_HD/R_code/my_function/MyPlot_v07.RData")



##################################################################
# fig3 B compartment overlap(有多少NAIR在B中)
##################################################################
hela_NAT_table.fix = hela_NAT_table
# hela_NAT_table.fix = hela_NAT_table[(hela_NAT_table$V3-hela_NAT_table$V2) >= 1000e3,]

NAT_in_B_ratio = rep(0,23)

for(chrom_index in c(1:23)){
  print(chrom_index)
  chrom_name = chrom.name(chrom_index)
  chrom_NAT_table = hela_NAT_table.fix[hela_NAT_table.fix$V1==chrom_name,]
  chrom_NAT.overlap = overlap(chrom_NAT_table,hela_AB_table,level_vetor = chrom_name)
  chrom_NAT.overlap.B = chrom_NAT.overlap[chrom_NAT.overlap$name.2=="B",]
  
  chrom_B_table = hela_AB_table[hela_AB_table$chrom_name==chrom_name & hela_AB_table$info=="B",]
  
  NAT_in_B_ratio[chrom_index] = sum(chrom_NAT.overlap.B$end - chrom_NAT.overlap.B$start) / sum(chrom_NAT.overlap[,3] - chrom_NAT.overlap[,2])
}

NAT_in_B_ratio.table = data.frame(chrom_name = c(paste("chr",c(1:22),sep = ""),"chrX"),
                                  NAT_in_B_ratio = NAT_in_B_ratio)

HG19_LEN=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
           135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,
           51304566,155270560,59373566,42999)

## 计算genome的平均值
# length_vector = HG19_LEN[1:23] / sum(HG19_LEN[1:23])
# NAT_in_B_ratio.mean = sum(NAT_in_B_ratio.table$NAT_in_B_ratio * length_vector)
NAT_in_B_ratio.mean = length(which(hela_AB_table$info=="B")) / nrow(hela_AB_table)

# barplot 
par(mar=c(3,3,2,1),family="Arial",font=1)
barplot(NAT_in_B_ratio.table$NAT_in_B_ratio,col = ifelse(NAT_in_B_ratio.table$NAT_in_B_ratio>=NAT_in_B_ratio.mean,yes = "#CC5C55",no = "#42B256"),border = F,xaxt="n",yaxt="n",ylim = c(0,1))
abline(h = NAT_in_B_ratio.mean,lwd=3,pch=10,lty=2,col="black")
axis(side = 2,at=seq(0,1,0.25),labels = F,cex.axis=3,lwd=3,tck=-0.05)



##################################################################
# fig3 NAIR能解释多少B
##################################################################
# 将染色体分成2类
# 1类是与核仁相近的，能够用NAIR解释B-compartment；
# 2类是与染色体相距较远的，NAIR不能够很好解释；

hela_NAT_table.fix = hela_NAT_table
# hela_NAT_table.fix = hela_NAT_table[(hela_NAT_table$V3-hela_NAT_table$V2) >= 1000e3,]

NAT_B_ratio = rep(0,23)
for(chrom_index in c(1:23)){
  print(chrom_index)
  chrom_name = chrom.name(chrom_index)
  chrom_NAT_table = hela_NAT_table.fix[hela_NAT_table.fix$V1==chrom_name,]
  chrom_NAT.overlap = overlap(chrom_NAT_table,hela_AB_table,level_vetor = chrom_name)
  chrom_NAT.overlap.B = chrom_NAT.overlap[chrom_NAT.overlap$name.2=="B",]
  
  chrom_B_table = hela_AB_table[hela_AB_table$chrom_name==chrom_name & hela_AB_table$info=="B",]
  
  NAT_B_ratio[chrom_index] = sum(chrom_NAT.overlap.B$end - chrom_NAT.overlap.B$start) / sum(chrom_B_table[,3] - chrom_B_table[,2])
}

NAT_B_ratio.table = data.frame(chrom_name = c(paste("chr",c(1:22),sep = ""),"chrX"),
                               NAT_B_ratio = NAT_B_ratio)

## 按照染色体从小到大排序
HG19_LEN=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
           135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,
           51304566,155270560,59373566,42999)

## 计算genome的平均值
length_vector = HG19_LEN[1:23] / sum(HG19_LEN[1:23])
NAT_B_ratio.mean = sum(NAT_B_ratio * length_vector)

# 按长度排序从小到大
# NAT_B_ratio.table$chrom_name = factor(x = NAT_B_ratio.table$chrom_name,levels = chrom.name.factor(order(HG19_LEN[1:23])))

# 按序号排序 1-23
NAT_B_ratio.table$chrom_name = factor(x = NAT_B_ratio.table$chrom_name,levels = chrom.name.factor(1:23))

# barplot 
par(mar=c(6,6,2,1),family="Arial",font=1)
barplot(NAT_B_ratio.table$NAT_B_ratio,col = ifelse(NAT_B_ratio.table$NAT_B_ratio>=NAT_B_ratio.mean,yes = "#FF6B00",no = "#7061CC"),border = F,xaxt="n",yaxt="n",ylim=c(0,0.8))
abline(h = NAT_B_ratio.mean,lwd=3,pch=10,lty=2,col="black")
axis(side = 2,at=seq(0,0.8,0.2),labels = F,cex.axis=3,lwd=3,tck=-0.05)


##################################################################
# fig3 signal strength
##################################################################
hela_NAT_table.fix = hela_NAT_table
hela_NAT_table.fix = cbind(hela_NAT_table.fix,rep(0,nrow(hela_NAT_table.fix)))
colnames(hela_NAT_table.fix) = c("chrom_name","start","end","peak_index","value","ratio")

# for(row.index in c(1:nrow(hela_NAT_table))){

for(row.index in c(1:nrow(hela_NAT_table))){  
  print(row.index)
  table.row = hela_NAT_table[row.index,]
  chrom_name = as.character(table.row[1,1])
  region_start = as.integer(table.row[1,2])
  region_end = as.integer(table.row[1,3])
  
  select_vector = hic_colsum_table$chrom_name == as.character(chrom_name) & hic_colsum_table$region_start >= region_start & hic_colsum_table$region_end <= region_end
  n_hic_colsum_vector = hic_colsum_table$n_hic_colsum[select_vector]
  g_hic_colsum_vector = hic_colsum_table$g_hic_colsum[select_vector]
  
  ratio = (sum(n_hic_colsum_vector) / hic_count.hela_n_hic.total) / (sum(g_hic_colsum_vector) / hic_count.hela_g_hic.total)
  
  hela_NAT_table.fix$ratio[row.index] = ratio
}

chrom_NAT_median = rep(1,23)
for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_NAT_median[chrom_index] = median(hela_NAT_table.fix$ratio[hela_NAT_table.fix$chrom_name==chrom_name])
}

par(mar=c(6,6,2,1),family="Arial",font=1)
barplot(chrom_NAT_median,col = ifelse(chrom_NAT_median>=median(hela_NAT_table.fix$ratio),yes = "#FF6B00",no = "#7061CC"),border = F,xaxt="n",yaxt="n",ylim=c(0,3))
abline(h = NAT_B_ratio.mean,lwd=3,pch=10,lty=2,col="black")
axis(side = 2,at=seq(0,0.8,0.2),labels = F,cex.axis=3,lwd=3,tck=-0.05)


ratio_vector = as.vector(NAT_B_ratio.table$NAT_B_ratio)

plot(x=ratio_vector,
     y=chrom_NAT_median,
     xlim=c(0,0.8),
     ylim=c(1.5,3),
     col=ifelse(ratio_vector >= NAT_B_ratio.mean,yes = "#FF6B00",no = "#7061CC"),
     pch=19,frame.plot = F,xaxt="n",yaxt="n",xlab="",ylab="",cex=2)

axis(side = 1,at=seq(0,0.8,0.2),labels = F,cex.axis=3,lwd=3,tck=-0.05)
axis(side = 2,at=seq(1.5,3,0.5),labels = F,cex.axis=3,lwd=3,tck=-0.05)

text(x=ratio_vector[ratio_vector>=NAT_B_ratio.mean],y=chrom_NAT_median[ratio_vector>=NAT_B_ratio.mean]+ 0.05,labels = NAT_B_ratio.table$chrom_name[ratio_vector>=NAT_B_ratio.mean],cex=2)






