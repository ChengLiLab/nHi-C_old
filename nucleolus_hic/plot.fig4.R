#############################################################################################################
# 2016-09-21 Howard MENG
# plot fig3 figures
# 1.论证NAIR 和TAD的统一性
# 2.NAIR的边界上有CTCF的结合peak
# 3.NAIR的
#############################################################################################################

###########################################################################
# fig3 NAIR-TAD plot
###########################################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")

hela_NAT = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
hela_NAD = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep = "\t")

hela_g_hic_insulation = read.table("~/menghw_HD/data_table/TAD_region/hela-g-hic_insulation_score_w1e6.bed",header = F,sep = "\t")
colnames(hela_g_hic_insulation) = c("chrom_name","start","end","score")

hela_n_hic_insulation = read.table("~/menghw_HD/data_table/TAD_region/hela-n-hic_insulation_score_w1e6.bed",header = F,sep = "\t")
colnames(hela_n_hic_insulation) = c("chrom_name","start","end","score")

hela_g_hic_TAD = read.table("~/menghw_HD/data_table/TAD_region/hela-g-hic_TAD.bed",header = F,sep = "\t")
colnames(hela_g_hic_TAD) = c("chrom_name","start","end")
hela_g_hic_TAD = hela_g_hic_TAD[(hela_g_hic_TAD$end - hela_g_hic_TAD$start) > 200e3,]

hela_n_hic_TAD = read.table("~/menghw_HD/data_table/TAD_region/hela-n-hic_TAD.bed",header = F,sep = "\t")
colnames(hela_n_hic_TAD) = c("chrom_name","start","end")
hela_n_hic_TAD = hela_n_hic_TAD[(hela_n_hic_TAD$end - hela_n_hic_TAD$start) > 0,]


chrom_index = 15
binsize = 20e3
load(file=sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_%d_%d.RData",chrom_index,binsize))

hela_NAT.fix = hela_NAT[hela_NAT$V3-hela_NAT$V2 >= 500e3,]

chrom_start = 52e6
chrom_end = 63e6
chrom_name = chrom.name(chrom_index)

g_hic_matrix.part = matrix.part(g_hic_matrix.norm,chrom_index,chrom_start,chrom_end,binsize)
n_hic_matrix.part = matrix.part(n_hic_matrix,chrom_index,chrom_start,chrom_end,binsize)

fig.facet <- layout(matrix(c(1:4),nrow = 4,byrow = FALSE),width = rep(10,4),heights = c(7,1,7,2))
layout.show(fig.facet)

par(mar=c(0.5,2,0.5,2))
plot.TAD(g_hic_matrix.part,ylim=c(0,100),maxBound = 0.97)

par(mar=c(0.5,2,0.5,2))
# plot.peak(df.peak = hela_NAT,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
plot.peak(df.peak = hela_NAT.fix,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

# par(mar=c(0.5,2,0.5,2))
# plot.peak(df.peak = hela_g_hic_TAD.overlap.raw.fix,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF0000",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

par(mar=c(0.5,2,0.5,2))
plot.TAD(n_hic_matrix.part,ylim=c(-100,0),maxBound = 0.958,mat.upper = F)

# par(mar=c(0.5,2,0.5,2))
# plot.barplot(hela_CTCF.fix,chrom_index,chrom_start,chrom_end)


fig.facet <- layout(matrix(c(1:5),nrow = 5,byrow = FALSE),width = rep(10,5),heights = c(3,2,3,2,1))
layout.show(fig.facet)

# hela-g-hic
par(mar=c(0.5,2,0.5,2))
plot.TAD(g_hic_matrix.part,ylim=c(0,50),maxBound = 0.97)

# hela-g-hic insulation score and TAD boundary
par(mar=c(0.5,2,0,2))
g_hic.score.value.part = hela_g_hic_insulation$score[hela_g_hic_insulation$chrom_name==chrom_name & hela_g_hic_insulation$start >= chrom_start & hela_g_hic_insulation$start <= chrom_end]
g_hic.score.start.part = hela_g_hic_insulation$start[hela_g_hic_insulation$chrom_name==chrom_name & hela_g_hic_insulation$start >= chrom_start & hela_g_hic_insulation$start <= chrom_end]
plot(g_hic.score.start.part,g_hic.score.value.part,type="l",frame.plot=F,xaxt="n",yaxt="n",lwd=3,col="#00685A",xlim=c(chrom_start,chrom_end),ylim=c(-1,1),xlab="",ylab="")
axis(side=2,at=c(-1,0,1),labels = F,lwd=3,tck=-0.2)

par(new=T)
g_hic.TAD.start.part = hela_g_hic_TAD$start[hela_g_hic_TAD$chrom_name==chrom_name & hela_g_hic_TAD$start >= chrom_start & hela_g_hic_TAD$start <= chrom_end]
g_hic.TAD.start.part = g_hic.TAD.start.part - binsize
plot(x=g_hic.TAD.start.part,y=rep(1,length(g_hic.TAD.start.part)),type="h",frame.plot=F,xaxt="n",yaxt="n",lwd=3,col="#1E0B5D",xlim=c(chrom_start,chrom_end),ylim=c(-2,1),xlab="",ylab="")

# hela-n-hic
par(mar=c(0.5,2,0.5,2))
plot.TAD(n_hic_matrix.part,ylim=c(0,50),maxBound = 0.958,mat.upper = T)

par(mar=c(0.5,2,0,2))
n_hic.score.value.part = hela_n_hic_insulation$score[hela_n_hic_insulation$chrom_name==chrom_name & hela_n_hic_insulation$start >= chrom_start & hela_n_hic_insulation$start <= chrom_end]
n_hic.score.value.part = n_hic.score.value.part / 3

n_hic.score.start.part = hela_n_hic_insulation$start[hela_n_hic_insulation$chrom_name==chrom_name & hela_n_hic_insulation$start >= chrom_start & hela_n_hic_insulation$start <= chrom_end]
plot(n_hic.score.start.part,n_hic.score.value.part,type="l",frame.plot=F,xaxt="n",yaxt="n",lwd=3,col="#00685A",xlim=c(chrom_start,chrom_end),ylim=c(-1,1),xlab="",ylab="")
axis(side=2,at=c(-1,0,1),labels = F,lwd=3,tck=-0.2)

par(new=T)
g_hic.TAD.start.part = hela_g_hic_TAD$start[hela_g_hic_TAD$chrom_name==chrom_name & hela_g_hic_TAD$start >= chrom_start & hela_g_hic_TAD$start <= chrom_end]
g_hic.TAD.start.part = g_hic.TAD.start.part - binsize
plot(x=g_hic.TAD.start.part,y=rep(1,length(g_hic.TAD.start.part)),type="h",frame.plot=F,xaxt="n",yaxt="n",lwd=3,col="#1E0B5D",xlim=c(chrom_start,chrom_end),ylim=c(-2,1),xlab="",ylab="")


# par(new=T)
# n_hic.TAD.start.part = hela_n_hic_TAD$start[hela_n_hic_TAD$chrom_name==chrom_name & hela_n_hic_TAD$start >= chrom_start & hela_n_hic_TAD$start <= chrom_end]
# n_hic.TAD.start.part = n_hic.TAD.start.part - binsize
# plot(x=n_hic.TAD.start.part,y=rep(1,length(n_hic.TAD.start.part)),type="h",frame.plot=F,xaxt="n",yaxt="n",lwd=3,col="#1E0B5D",xlim=c(chrom_start,chrom_end),ylim=c(-2,1),xlab="",ylab="")


par(mar=c(1,2,0,2))
plot.peak(df.peak = hela_NAT.fix,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

dev.off()




###########################################################################
# fig3 TAD length distribution
###########################################################################
rm(list=ls())
hela_g_hic_TAD = read.table("~/menghw_HD/data_table/TAD_region/hela-g-hic_TAD.bed",header = F,sep = "\t")
colnames(hela_g_hic_TAD) = c("chrom_name","start","end")

hela_g_hic_TAD$start = as.numeric(hela_g_hic_TAD$start)
hela_g_hic_TAD$end = as.numeric(hela_g_hic_TAD$end)

hela_g_hic_TAD.fix = hela_g_hic_TAD[(hela_g_hic_TAD$end-hela_g_hic_TAD$start)>0,]
hela_g_hic_TAD.fix = hela_g_hic_TAD[(hela_g_hic_TAD$end-hela_g_hic_TAD$start)>= 120e3,]


sum(hela_g_hic_TAD.fix$end - hela_g_hic_TAD.fix$start)
quantile((hela_g_hic_TAD.fix$end - hela_g_hic_TAD.fix$start))
hist(log10(hela_g_hic_TAD.fix$end - hela_g_hic_TAD.fix$start))

dev.off()

TAD_length = (hela_g_hic_TAD.fix$end - hela_g_hic_TAD.fix$start)
TAD_median = median((hela_g_hic_TAD.fix$end - hela_g_hic_TAD.fix$start))


par(mar=c(2,2,2,2),family="Arial",font=1)
hist(log10(TAD_length),breaks = seq(4.5,8,0.1),xlim = c(4.5,7),ylim=c(0,600),xlab = "",ylab = "",main="",cex.axis=2,cex.lab=2,lwd=3,col="#00685A",xaxt="n",yaxt="n")

axis(side=1,at=seq(4.5,7,0.5),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
axis(side=2,at=seq(0,600,150),labels = F,cex.axis=3,lwd=3,tck=-0.05)

# text(x=4.8,y=seq(0,1200,400),labels = seq(0,1200,400),cex = 3,xpd=T,pos = 2)
# text(x=3.8,y=600,labels = "Count",cex = 3,xpd=T,srt=90)

abline(v=log10(TAD_median),lwd=5,pch=10,lty=2,col="black")
# text(x=6.5,y=1100,sprintf("median=%dKb",ceiling(TAD_median %/% 1000)),pos=1,cex=2.5)
# box(lwd=3)


###########################################################################
# fig3 NAIR length distribution
###########################################################################
rm(list=ls())
hela_g_hic_NAIR = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
colnames(hela_g_hic_NAIR) = c("chrom_name","start","end","index","info")

hela_g_hic_NAIR$start = as.numeric(hela_g_hic_NAIR$start)
hela_g_hic_NAIR$end = as.numeric(hela_g_hic_NAIR$end)

hela_g_hic_NAIR.fix = hela_g_hic_NAIR[(hela_g_hic_NAIR$end-hela_g_hic_NAIR$start)>0,]
hela_g_hic_NAIR.fix = hela_g_hic_NAIR[(hela_g_hic_NAIR$end-hela_g_hic_NAIR$start)>= 120e3,]


sum(hela_g_hic_NAIR.fix$end - hela_g_hic_NAIR.fix$start)
quantile((hela_g_hic_NAIR.fix$end - hela_g_hic_NAIR.fix$start))
hist(log10(hela_g_hic_NAIR.fix$end - hela_g_hic_NAIR.fix$start))

dev.off()

NAIR_length = (hela_g_hic_NAIR.fix$end - hela_g_hic_NAIR.fix$start)
NAIR_median = median((hela_g_hic_NAIR.fix$end - hela_g_hic_NAIR.fix$start))

par(mar=c(2,2,2,2),family="Arial",font=1)
hist(log10(NAIR_length),breaks = seq(4.5,8,0.1),xlim = c(4.5,7.5),ylim=c(0,60),xlab = "",ylab = "",main="",cex.axis=2,cex.lab=2,lwd=3,col="#FF7304",xaxt="n",yaxt="n")

axis(side=1,at=seq(4.5,7.5,0.5),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
axis(side=2,at=seq(0,60,15),labels = F,cex.axis=3,lwd=3,tck=-0.05)
abline(v=log10(NAIR_median),lwd=5,pch=10,lty=2,col="black")


###########################################################################
# fig3 barplot overlap between TAD & NAIR
###########################################################################
par(mar=c(3,3,2,1),family="Arial",font=1)
barplot(c(222,149),col = c("#00685A","#FF7304"),border = F,xaxt="n",yaxt="n",ylim = c(0,300))
axis(side = 2,at=seq(0,300,100),labels = F,cex.axis=3,lwd=5,tck=-0.1)
dev.off()


###########################################################################
# fig3 TAD, NAIR boundary -> insulation score
###########################################################################
#---- 读取数据----------------------------


TAD_boundary <- read.table('~/menghw_HD/data_table/TAD_region/hela-g-hic_TAD.bed', header = T, as.is = T)
insulationScore <- read.table("~/menghw_HD/data_table/TAD_region/hela-g-hic_insulation_score_w1e6.bed", header = F, as.is = T)
NAIR <- read.table('~/menghw_HD/data_table/Hela_NAT.bed', header = F, as.is = T)

# turn chrX into chr23:
head(insulationScore) # chr1-23
colnames(insulationScore) <- c('chr', 'start', 'end', 'insulation_score')
NAIR[, 4] <- NAIR[, 3] - NAIR[, 2]
colnames(NAIR) <- c('chr', 'start', 'end', 'length')
NAIR[ NAIR$chr == 'chrX', 1 ] = 'chr23'
table(NAIR$chr)
#NAIR.boundary <- rbind(as.matrix(NAIR[, c(1, 2)]), as.matrix(NAIR[, c(1, 3)]))

#--- Function to generate boundary insulation score matrix from a given n x 2 data.frame----------------
# input data.frame format:  chr_name location
#                     e.g.    chr1     20000
generate_boundary_insulatin_score_matrix <- function( dataframe, insulationScore,  window = '2mb'  ){
  
  # 生成一个矩阵，每一行都是一个boundary,包括了其上下游个1Mb的insulationScore的信息
  boundary_insulationScore <- matrix(nrow = dim(dataframe)[1],  ncol = 101 )
  
  # 每一个循环处理一个TAD boundary
  # 第一步：取出insulationScore
  # 第二步：放入新的matrix中
  # 分三种情况：boundary在N, C末端，或者在中间位置。
  for(i in seq(1, dim(dataframe)[1])){
    chr <- dataframe[i, 1] # chromosome name
    sub_insulationScore <- insulationScore[insulationScore$chr == chr, ] # 此boundary所在的染色体的score矩阵
    index = which( dataframe[i, 2]  <= as.numeric(sub_insulationScore$end) )[1]  # 此boundary所在的bin的位置
    if(index < 51){ # N terminal 
      score <- sub_insulationScore[ ( 1 :  (index + 50)) , ]$insulation_score
      boundary_insulationScore[i, ( 51 - index + 1) : 101 ] <- score
      
    }else if( (dim(sub_insulationScore)[1] - index) < 50){ # C terminal
      score <- sub_insulationScore[ ( (index - 50) : dim(sub_insulationScore)[1] ) , ]$insulation_score 
      boundary_insulationScore[i, 1: ( 51 + (dim(sub_insulationScore)[1] - index) )] <- score
      
    }else{ # In the middle
      boundary_insulationScore[i,] <- sub_insulationScore[ ((index - 50) : (index + 50)), ]$insulation_score  
    }
  }
  return(boundary_insulationScore)
}

#----simulation for random conditions-------
#chr_length <- read.table('~/R/data/chr.length.txt', header = T, as.is = T)
#chr_length[23, 1] <- 'chr23'
#chr_length <- chr_length[1:23, ]
#chr_length[23, 1] = 'chr23'
#chr_num_NAIR <- table(NAIR.boundary$chr) #每条染色体的NAT boundary数目

#生成随机的NAT boundary信息：
#random_NAIR_boundary <- matrix(nrow = dim(NAIR.boundary)[1], ncol = 2)

#start_index = 1
#for(i in seq(1, 23)){
# chr_name <- chr_length[i, 1]
# chr_len <- chr_length[i, 2]
# sampling_num <- sample(chr_len, chr_num_NAIR[chr_name])
# random_NAIR_boundary[start_index : (start_index + length(sampling_num ) - 1) , 1 ] = rep(chr_name, length(sampling_num))
# random_NAIR_boundary[start_index : (start_index + length(sampling_num) - 1),  2 ] = sampling_num
# start_index = start_index + length(sampling_num) 
# }

#colnames(random_NAIR_boundary) <- c('chr', 'start')

#write.table(random_NAIR_boundary, file = '~/R/data/nucleulos/random_NAIR_boundary.txt', quote = F, row.names = F)
#write.table(NAIR.boundary, file = '~/R/data/nucleulos/NAIR_boundary.txt', quote = F, row.names = F)
random_NAIR_boundary <- read.table('~/menghw_HD/data_table/boundary_plot/random_NAIR_boundary.txt', header = T, as.is = T)
NAIR.boundary <- read.table('~/menghw_HD/data_table/boundary_plot/NAIR_boundary.txt', header = T, as.is = T)

#---- boundary insulationscores-------------
TAD_boundary_insulationScore <- generate_boundary_insulatin_score_matrix(TAD_boundary, insulationScore)
NAIR_boundary_insulationScore <- generate_boundary_insulatin_score_matrix(NAIR.boundary , insulationScore)
random_NAT_boundary_insulationScore <- generate_boundary_insulatin_score_matrix(random_NAIR_boundary, insulationScore)

#--plot --------------------------------

score = function(list){ return (median(as.numeric(list), na.rm = T)) }

#pdf('~/R/data/nucleulos/figures/boundary_insulation_score2.pdf')


layout(matrix(1))
par(mar=c(2,2,2,1))
plot.new()
plot.window(xlim = c(0, 110), ylim =c(-0.45, 0.3)   )
lines(1:101, apply(random_NAT_boundary_insulationScore, 2, score), col = '#003333', lwd = 2.5)
lines(1:101, apply(TAD_boundary_insulationScore, 2, score), col = '#00685A', lwd = 5)
lines(1:101, apply(NAIR_boundary_insulationScore, 2, score), col = '#FF7304', lwd = 5)

# axis(side = 2, lwd = 2, at = c(-0.4, -0.2, 0, 0.2), labels = F, tck= -0.01)

axis(side = 1, lwd = 3,at = c(0,25,50,75,100), tck= -0.05,labels = F)
axis(side = 2, lwd = 3, at = c(-0.4, -0.2, 0, 0.2), tck= -0.05,labels = F)

###########################################################################
# fig3 TAD, NAIR boundary -> histone modification
###########################################################################
rm(list=ls())
# plot the distribution of Chip seq data at TAD(NAIR) boundaries
#-- import data ------------
NAIR.boundary <- read.table('~/menghw_HD/data_table/boundary_plot/NAIR_boundary.txt', header = T, as.is = T)
TAD_boundary <- read.table("~/menghw_HD/data_table/TAD_region/hela-g-hic_TAD.bed", header = F, as.is = T)
random_NAIR_boundary <- read.table("~/menghw_HD/data_table/boundary_plot/random_NAIR_boundary.txt", header = T, as.is = T)

Chip.CTCF <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_CTCF_ENCFF000BAJ_coverage_10000.bed', header = F, as.is = T)
colnames(Chip.CTCF) <- c('chr', 'start', 'end', 'reads_number', 'coverage', 'length', 'ratio')

Chip.H3K27me3 <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_H3K27me3_ENCFF000BBS_coverage_10000.bed', header = F, as.is = T)
colnames(Chip.H3K27me3) <- c('chr', 'start', 'end', 'reads_number', 'coverage', 'length', 'ratio')

Chip.H3K36me3 <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_H3K36me3_ENCFF000BCA_coverage_10000.bed', header = F, as.is = T)
colnames(Chip.H3K36me3) <- c('chr', 'start', 'end', 'reads_number', 'coverage', 'length', 'ratio')

Chip.H3K4me3 <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_H3K4me3_ENCFF000BCO_coverage_10000.bed', header = F, as.is = T)
colnames(Chip.H3K4me3) <- c('chr', 'start', 'end', 'reads_number', 'coverage', 'length', 'ratio')

Chip.H3K9me3 <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_H3K9me3_ENCFF000BBG_coverage_10000.bed', header = F, as.is = T)
colnames(Chip.H3K9me3) <- c('chr', 'start', 'end', 'reads_number', 'coverage', 'length', 'ratio')

Chip.H4K20me1 <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_H4K20me1_ENCFF000BDC_coverage_10000.bed', header = F, as.is = T)
colnames(Chip.H4K20me1) <- c('chr', 'start', 'end', 'reads_number', 'coverage', 'length', 'ratio')

# chrX --> chr23, discard chrY -------
Change.chrXY <- function(dataframe){
  dataframe[dataframe$chr == 'chrX', 1] = 'chr23'
  dataframe = dataframe[ dataframe$chr != 'chrY', ]
  return(dataframe)
}

Chip.CTCF <- Change.chrXY(Chip.CTCF)
Chip.H4K20me1 <- Change.chrXY(Chip.H4K20me1)
Chip.H3K9me3 <- Change.chrXY(Chip.H3K9me3)
Chip.H3K4me3 <- Change.chrXY(Chip.H3K4me3)
Chip.H3K36me3 <- Change.chrXY(Chip.H3K36me3)
Chip.H3K27me3 <- Change.chrXY(Chip.H3K27me3)

# functions 
generate_boundary_chip_readsnum_matrix <- function( dataframe, ChipSeq,  window = '1mb'  ){
  
  # 生成一个矩阵，每一行都是一个boundary,包括了其上下游个1Mb的ChipSeq的信息
  boundary_ChipSeq <- matrix(nrow = dim(dataframe)[1],  ncol = 101 )
  
  # 每一个循环处理一个TAD boundary
  # 第一步：取出ChipSeq
  # 第二步：放入新的matrix中
  # 分三种情况：boundary在N, C末端，或者在中间位置。
  for(i in seq(1, dim(dataframe)[1])){
    chr <- dataframe[i, 1] # chromosome name
    sub_ChipSeq <- ChipSeq[ChipSeq$chr == chr, ] # 此boundary所在的染色体的score矩阵
    index = which( dataframe[i, 2]  <= as.numeric(sub_ChipSeq$end) )[1]  # 此boundary所在的bin的位置
    if(index < 51){ # N terminal 
      score <- sub_ChipSeq[ ( 1 :  (index + 50)) , ]$reads_number
      boundary_ChipSeq[i, ( 51 - index + 1) : 101 ] <- score
      
    }else if( (dim(sub_ChipSeq)[1] - index) < 50){ # C terminal
      score <- sub_ChipSeq[ ( (index - 50) : dim(sub_ChipSeq)[1] ) , ]$reads_number 
      boundary_ChipSeq[i, 1: ( 51 + (dim(sub_ChipSeq)[1] - index) )] <- score
      
    }else{ # In the middle
      boundary_ChipSeq[i,] <- sub_ChipSeq[ ((index - 50) : (index + 50)), ]$reads_number  
    }
  }
  return(boundary_ChipSeq)
}

# NAIR.boundary calculatioin----
random.NAIR.boundary.Chip.CTCF <- generate_boundary_chip_readsnum_matrix(random_NAIR_boundary, Chip.CTCF)
NAIR.boundary.Chip.CTCF <- generate_boundary_chip_readsnum_matrix(NAIR.boundary, Chip.CTCF)

random.NAIR.boundary.Chip.H3K27me3 <- generate_boundary_chip_readsnum_matrix(random_NAIR_boundary, Chip.H3K27me3)
NAIR.boundary.Chip.H3K27me3 <- generate_boundary_chip_readsnum_matrix(NAIR.boundary, Chip.H3K27me3)

random.NAIR.boundary.Chip.H3K36me3 <- generate_boundary_chip_readsnum_matrix(random_NAIR_boundary, Chip.H3K36me3)
NAIR.boundary.Chip.H3K36me3 <- generate_boundary_chip_readsnum_matrix(NAIR.boundary, Chip.H3K36me3)

random.NAIR.boundary.Chip.H3K4me3 <- generate_boundary_chip_readsnum_matrix(random_NAIR_boundary, Chip.H3K4me3)
NAIR.boundary.Chip.H3K4me3 <- generate_boundary_chip_readsnum_matrix(NAIR.boundary, Chip.H3K4me3)

random.NAIR.boundary.Chip.H3K9me3 <- generate_boundary_chip_readsnum_matrix(random_NAIR_boundary, Chip.H3K9me3)
NAIR.boundary.Chip.H3K9me3 <- generate_boundary_chip_readsnum_matrix(NAIR.boundary, Chip.H3K9me3)

random.NAIR.boundary.Chip.H4K20me1 <- generate_boundary_chip_readsnum_matrix(random_NAIR_boundary, Chip.H4K20me1)
NAIR.boundary.Chip.H4K20me1 <- generate_boundary_chip_readsnum_matrix(NAIR.boundary, Chip.H4K20me1)

# TAD.boundary calculation ------
TAD_boundary_Chip.CTCF <- generate_boundary_chip_readsnum_matrix(TAD_boundary, Chip.CTCF)
TAD_boundary_Chip.H3K4me3 <- generate_boundary_chip_readsnum_matrix(TAD_boundary, Chip.H3K4me3)
TAD_boundary_Chip.H3K9me3 <- generate_boundary_chip_readsnum_matrix(TAD_boundary, Chip.H3K9me3)
TAD_boundary_Chip.H3K27me3 <- generate_boundary_chip_readsnum_matrix(TAD_boundary, Chip.H3K27me3)
TAD_boundary_Chip.H3K36me3 <- generate_boundary_chip_readsnum_matrix(TAD_boundary, Chip.H3K36me3)
TAD_boundary_Chip.H4K20me1 <- generate_boundary_chip_readsnum_matrix(TAD_boundary, Chip.H4K20me1)

# Plot ---------------
score = function(list){ return (mean(as.numeric(list), na.rm = T)) } # mean than median

plot.chip <- function(random.data, sample.data, range = c(50, 120), lwd  = 5, at.y = c(50, 75, 100),col1 = '#003333', col2= '#2219b2', main = 'main', tck = -0.05){
  plot.new()
  plot.window(xlim = c(0, 110), ylim = range  )
  # title(main = main)
  lines(1:101, apply(random.data, 2, score), col = col1, lwd = lwd * 0.5)
  lines(1:101, apply(sample.data, 2, score), col = col2, lwd = lwd)
  axis(side = 1, lwd = 3,at = c(0,25,50,75,100), tck= tck,labels = F)
  axis(side = 2, lwd = 3, at = at.y, tck= tck,labels = F)
}

# NAIR plot 1 -------
#pdf('~/R/data/nucleulos/figures/NAIR_boundary_chip.pdf')
#png('~/R/data/nucleulos/figures/NAIR_boundary_chip.png')
# layout(matrix(c(1,2,3,4, 5, 6), 2,3))
# par(oma = c(4 , 4, 4 ,4))
# par(mar = c(2, 0.5, 2, 0.5))

par(mar=c(2,2,0.5,0.5))
#### active
# H3K4me3
plot.chip(random.NAIR.boundary.Chip.H3K4me3, NAIR.boundary.Chip.H3K4me3,col2 = '#497B1F',range = c(50, 110),at.y = seq(50, 110, 20))
# H4K20me1 
plot.chip(random.NAIR.boundary.Chip.H4K20me1, NAIR.boundary.Chip.H4K20me1, main = 'H4K20me1', col2 = '#EA7B55',range = c(50, 130), at.y = seq(50, 130, 20))
# H3K36me3
plot.chip(random.NAIR.boundary.Chip.H3K36me3, NAIR.boundary.Chip.H3K36me3, range = c(70, 150), main = 'H3K36me3', col2 = '#D14C74', at.y = seq(70, 150, 20))

#### repressive
# CTCF
plot.chip(random.NAIR.boundary.Chip.CTCF, NAIR.boundary.Chip.CTCF, range = c(60, 120),main = 'CTCF',col2 = "#8A4723",at.y = seq(60,120,20))
# H3K9me3
plot.chip(random.NAIR.boundary.Chip.H3K9me3, NAIR.boundary.Chip.H3K9me3, range = c(100, 350), main = 'H3K9me3', col2 = '#1F2D5D', at.y = seq(100,360,80))
# H3K27me3
plot.chip(random.NAIR.boundary.Chip.H3K27me3, NAIR.boundary.Chip.H3K27me3, range = c(100, 190), main = 'H3K27me3', col2 = '#997CAB', at.y = seq(100,190,30))

# #dev.off()
# 
# # NAIR plot 2 --------
# plot.chip2 <- function(random.data, sample.data, range = c(50, 120), lwd  = 3, at.y = c(50, 75, 100),
#                        col1 = '#FF9700', col2= '#2219b2', main = 'main', tck = -0.01){
#   plot.new()
#   plot.window(xlim = c(0, 110), ylim = range  )
#   lines(1:101, apply(random.data, 2, score), col = col1, lwd = lwd * 0.34)
#   lines(1:101, apply(sample.data, 2, score), col = col2, lwd = lwd * 0.67)
#   axis(side = 1, lwd = 2,at = c(0, 51, 100), tck= tck, labels = F)
#   axis(side = 2, lwd = 2, at = at.y, tck= tck, labels = F)
# }
# 
# # NAIR plot
# # pdf('~/R/data/nucleulos/figures/NAIR_boundary_chip2.pdf')
# #png('~/R/data/nucleulos/figures/NAIR_boundary_chip.png')
# layout(matrix(c(1,2,3,4, 5, 6), 2,3))
# par(oma = c(4 , 4, 4 ,4))
# par(mar = c(2, 0.5, 2, 0.5))
# 
# plot.chip2(random.NAIR.boundary.Chip.H3K4me3, NAIR.boundary.Chip.H3K4me3, main = 'H3K4me3', range = c(50, 110), col2 = '#FB3F51', at.y = c(50, 75, 100))
# plot.chip2(random.NAIR.boundary.Chip.H3K9me3, NAIR.boundary.Chip.H3K9me3, range = c(100, 350), main = 'H3K9me3', col2 = '#D235D2', at.y = c(100, 200, 300))
# plot.chip2(random.NAIR.boundary.Chip.H4K20me1, NAIR.boundary.Chip.H4K20me1, main = 'H4K20me1', col2 = '#36D695', at.y = c(50, 80, 110))
# plot.chip2(random.NAIR.boundary.Chip.H3K27me3, NAIR.boundary.Chip.H3K27me3, range = c(100, 190), main = 'H3K27me3', col2 = '#33CCCC', at.y = c(100, 140, 180))
# plot.chip2(random.NAIR.boundary.Chip.H3K36me3, NAIR.boundary.Chip.H3K36me3, range = c(70, 150), main = 'H3K36me3', col2 = '#AC3BD4', at.y = c(80, 110, 140))
# plot.chip2(random.NAIR.boundary.Chip.CTCF, NAIR.boundary.Chip.CTCF, range = c(60, 120),main = 'CTCF', at.y = c(60, 90, 120))
# 
# #dev.off()
# 
# 
# layout(matrix(1))
# plot.new()
# # TAD plot
# plot.chip(random.NAIR.boundary.Chip.CTCF, TAD_boundary_Chip.CTCF)
# plot.chip(random.NAIR.boundary.Chip.H3K4me3, TAD_boundary_Chip.H3K4me3)
# plot.chip(random.NAIR.boundary.Chip.H3K9me3, TAD_boundary_Chip.H3K9me3, range = c(100, 200))
# plot.chip(random.NAIR.boundary.Chip.H4K20me1, TAD_boundary_Chip.H4K20me1)
# plot.chip(random.NAIR.boundary.Chip.H3K27me3, TAD_boundary_Chip.H3K27me3, range = c(100, 200))
# plot.chip(random.NAIR.boundary.Chip.H3K36me3, TAD_boundary_Chip.H3K36me3, range = c(70, 150))
# 
# 
# 
# plot.window(xlim = c(0, 110), ylim =c(63, 116) )
# title(main = main)
# lines(1:101, apply(random_NAIR_boundary_Chip.CTCF, 2, score), col = '#FF9700', lwd = 5)
# lines(1:101, apply(NAIR.boundary_Chip.CTCF, 2, score), col = '#2219b2', lwd = 5)
# 
# axis(side = 1, lwd = 2, at = c(0, 51, 100))
# axis(side = 2, lwd = 2)


###########################################################################
# fig3 plot outside NAIR region and inside NAIR region
###########################################################################
rm(list=ls())

#########################################################
## function 
#########################################################
chrom.name <- function(chrom_index){
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }else if(chrom_index == 25){
    chrom_name = "rDNA"
  }
}

make.boundary.matrix <- function(NAIR_boundary_table,coverage_table,cover_cutoff=0.5,window_size=10e3,inside_length=500e3,outside_length=500e3){
  # NAIR_boundary_table BED 
  # coverage_table BED file from bedtools
  
  ################### initial the plot matrix
  inside_matrix.ncol = inside_length %/% window_size
  outside_matrix.ncol = outside_length %/% window_size
  inside_matrix.nrow = nrow(NAIR_boundary_table)
  outside_matrix.nrow = nrow(NAIR_boundary_table)
  
  inside_matrix = matrix(data=NA,nrow = inside_matrix.nrow,ncol = inside_matrix.ncol)
  outside_matrix = matrix(data=NA,nrow = outside_matrix.nrow,ncol = outside_matrix.ncol)
  
  total_index = 1
  
  for(chrom_index in c(1:23)){
    chrom_name = chrom.name(chrom_index)
    print(chrom_name)
    chrom.NAIR.table = NAIR_boundary_table[NAIR_boundary_table$chrom_name==chrom_name,]
    chrom.coverage.table = coverage_table[coverage_table$chrom_name==chrom_name,]
    
    if(nrow(chrom.NAIR.table) >= 1){
      for(chrom.NAIR.index in c(1:nrow(chrom.NAIR.table))){
        chrom.NAIR.row = chrom.NAIR.table[chrom.NAIR.index,]
        NAIR.length = sum(chrom.NAIR.row[,3] - chrom.NAIR.row[,2])
        NAIR.length.index = NAIR.length %/% window_size
        
        ################## calculate inside information
        if(NAIR.length.index > inside_matrix.ncol){
          NAIR.length.index = inside_matrix.ncol
        }
        
        NAIR.start.index = chrom.NAIR.row[,2] %/% window_size
        table.start.index = NAIR.start.index + 1
        table.end.index = table.start.index + NAIR.length.index - 1 
        inside_matrix[total_index,c(1:(table.end.index - table.start.index + 1))] = as.vector(chrom.coverage.table$reads_count[c(table.start.index:table.end.index)])
        
        ################## calculate outside information  
        if(chrom.NAIR.index == 1){
          NAIR.outside.length = chrom.NAIR.row[,2]
        }else{
          NAIR.outside.length = chrom.NAIR.table[chrom.NAIR.index,2] - chrom.NAIR.table[chrom.NAIR.index-1,3]
        }
        NAIR.outside.length.index = NAIR.outside.length %/% window_size
        if(NAIR.outside.length.index > 0){
          if(NAIR.outside.length.index > outside_matrix.ncol){
            NAIR.outside.length.index = outside_matrix.ncol
          }
          
          outside.table.end.index =  NAIR.start.index
          outside.table.start.index = outside.table.end.index - NAIR.outside.length.index + 1
          outside_matrix[total_index,c((50 - NAIR.outside.length.index + 1):50)] = as.vector(chrom.coverage.table$reads_count[c(outside.table.start.index:outside.table.end.index)])
        }
        total_index = total_index + 1
      } 
    }
  }
  outside.inside.matrix = cbind(outside_matrix,inside_matrix)
  return(outside.inside.matrix)
}


plot.NAIR.region <- function(plot.matrix,plot.ylim=c(100,120),method="median"){
  plot.matrix.nrow = dim(plot.matrix)[1]
  plot.matrix.ncol = dim(plot.matrix)[2]
  median.vector = rep(NA,plot.matrix.ncol)
  mean.vector = rep(NA,plot.matrix.ncol)
  
  for(col.index in c(1:plot.matrix.ncol)){
    col.vector = plot.matrix[!is.na(plot.matrix[,col.index]),col.index]
    median.vector[col.index] = median(col.vector)
    mean.vector[col.index] = mean(col.vector)
  }
  # plot.ylim = c(100,115)
  # plot(x=c(1:100),y=mean.vector,type='l')
  if(method == "median"){
    plot(x=c(1:100),y=median.vector,type='l',lwd=5,col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=plot.ylim,frame.plot=F)  
  }else if(method == "mean"){
    plot(x=c(1:100),y=mean.vector,type='l',lwd=5,col="red",xaxt="n",yaxt="n",xlab="",ylab="",ylim=plot.ylim,frame.plot=F)  
  }
  
  box(lwd=3)
  # abline(v=50,lwd=3,col="#FF7304")
}

# save.image(file="~/menghw_HD/R_code/my_function/fig4.plot.RData")
load(file="~/menghw_HD/R_code/my_function/fig4.plot.RData")

hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header =F,sep = '\t')
colnames(hela_NAIR_table) = c("chrom_name","start","end","peak_index","info")

hela_CTCF <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_CTCF_ENCFF000BAJ_coverage_10000.bed', header = F, as.is = T)
colnames(hela_CTCF) = c('chrom_name', 'start', 'end', 'reads_count', 'coverage', 'length', 'ratio')

hela_H3K9me3 <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_H3K9me3_ENCFF000BBG_coverage_10000.bed', header = F, as.is = T)
colnames(hela_H3K9me3) = c('chrom_name', 'start', 'end', 'reads_count', 'coverage', 'length', 'ratio')

hela_H3K4me3 <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_H3K4me3_ENCFF000BCO_coverage_10000.bed', header = F, as.is = T)
colnames(hela_H3K4me3) = c('chrom_name', 'start', 'end', 'reads_count', 'coverage', 'length', 'ratio')

hela_H3K36me3 <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_H3K36me3_ENCFF000BCA_coverage_10000.bed', header = F, as.is = T)
colnames(hela_H3K36me3) = c('chrom_name', 'start', 'end', 'reads_count', 'coverage', 'length', 'ratio')

hela_H4K20me1 <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_H4K20me1_ENCFF000BDC_coverage_10000.bed', header = F, as.is = T)
colnames(hela_H4K20me1) = c('chrom_name', 'start', 'end', 'reads_count', 'coverage', 'length', 'ratio')

hela_H3K27me3 <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_H3K27me3_ENCFF000BBS_coverage_10000.bed', header = F, as.is = T)
colnames(hela_H3K27me3) = c('chrom_name', 'start', 'end', 'reads_count', 'coverage', 'length', 'ratio')

hela_POL2 <- read.table('~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_POL2_ENCFF000PEL_coverage_10000.bed', header = F, as.is = T)
colnames(hela_POL2) = c('chrom_name', 'start', 'end', 'reads_count', 'coverage', 'length', 'ratio')


plot.matrix.CTCF = make.boundary.matrix(NAIR_boundary_table = hela_NAIR_table,coverage_table = hela_CTCF)
plot.matrix.H3K4me3 = make.boundary.matrix(NAIR_boundary_table = hela_NAIR_table,coverage_table = hela_H3K4me3)
plot.matrix.H3K9me3 = make.boundary.matrix(NAIR_boundary_table = hela_NAIR_table,coverage_table = hela_H3K9me3)
plot.matrix.H3K27me3 = make.boundary.matrix(NAIR_boundary_table = hela_NAIR_table,coverage_table = hela_H3K27me3)
plot.matrix.H3K36me3 = make.boundary.matrix(NAIR_boundary_table = hela_NAIR_table,coverage_table = hela_H3K36me3)
plot.matrix.H4K20me1 = make.boundary.matrix(NAIR_boundary_table = hela_NAIR_table,coverage_table = hela_H4K20me1)
plot.matrix.POL2 = make.boundary.matrix(NAIR_boundary_table = hela_NAIR_table,coverage_table = hela_POL2)



## Hela-CTCF
plot.NAIR.region(plot.matrix.CTCF,plot.ylim = c(50,200),method = "mean")
axis(side=1,at=seq(0,100,25),labels = F,lwd=3,tck=-0.05)
axis(side=2,at=seq(50,200,50),labels = F,lwd=3,tck=-0.05)
rect(50,45,110,205,col="#FF730455",border = F)

## Hela-H3K4me3
plot.NAIR.region(plot.matrix.H3K4me3,plot.ylim = c(35,50),method = "median")
axis(side=1,at=seq(0,100,25),labels = F,lwd=3,tck=-0.05)
axis(side=2,at=seq(35,50,5),labels = F,lwd=3,tck=-0.05)
rect(50,30,110,55,col="#FF730455",border = F)

## Hela-H3K9me3
plot.NAIR.region(plot.matrix.H3K9me3,plot.ylim = c(95,120),method = "median")
axis(side=1,at=seq(0,100,25),labels = F,lwd=3,tck=-0.05)
axis(side=2,at=seq(95,120,10),labels = F,lwd=3,tck=-0.05)
rect(50,90,110,125,col="#FF730455",border = F)

## Hela-H3K27me3
plot.NAIR.region(plot.matrix.H3K27me3,plot.ylim = c(70,140),method="median")
axis(side=1,at=seq(0,100,25),labels = F,lwd=3,tck=-0.05)
axis(side=2,at=seq(70,140,20),labels = F,lwd=3,tck=-0.05)
rect(50,65,110,145,col="#FF730455",border = F)

## Hela-H3K36me3
plot.NAIR.region(plot.matrix.H3K36me3,plot.ylim = c(50,100),method="median")
axis(side=1,at=seq(0,100,25),labels = F,lwd=3,tck=-0.05)
axis(side=2,at=seq(50,100,15),labels = F,lwd=3,tck=-0.05)
rect(50,45,110,105,col="#FF730455",border = F)

## Hela-H4K20me1 *
plot.NAIR.region(plot.matrix.H4K20me1,plot.ylim = c(45,60),method="median")


## CTCF site * 
CTCF_site_outside_matrix = as.matrix(read.table("~/menghw_HD/data_table/fig4/outside_CTCF_site.matrix",sep = ",",header = F))
CTCF_site_inside_matrix = as.matrix(read.table("~/menghw_HD/data_table/fig4/inside_CTCF_site.matrix",sep = ",",header = F))
plot.matrix.CTCF_site = cbind(CTCF_site_outside_matrix,CTCF_site_inside_matrix)
plot.NAIR.region(plot.matrix.CTCF_site,plot.ylim = c(15,20),method = "median")


## POL2
plot.NAIR.region(plot.matrix.POL2,plot.ylim = c(25,55),method="median")
axis(side=1,at=seq(0,100,25),labels = F,lwd=3,tck=-0.05)
axis(side=2,at=seq(25,55,10),labels = F,lwd=3,tck=-0.05)
rect(50,20,110,60,col="#FF730455",border = F)



