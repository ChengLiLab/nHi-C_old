# plot the distribution of Chip seq data at TAD(NAIR) boundaries


#-- import data ------------
NAIR.boundary <- read.table('~/R/data/nucleulos/NAIR_boundary.txt', header = T, as.is = T)
TAD_boundary <- read.table("~/R/data/nucleulos/clean_data/boundary_hela_1m_cleanup.bed", header = F, as.is = T)
random_NAIR_boundary <- read.table('~/R/data/nucleulos/random_NAIR_boundary.txt', header = T, as.is = T)

Chip.CTCF <- read.table('~/R/data/chipseq/hela/Hela-ChIP_CTCF_ENCFF000BAJ_coverage_10000.bed', header = F, as.is = T)
colnames(Chip.CTCF) <- c('chr', 'start', 'end', 'reads_number', 'coverage', 'length', 'ratio')

Chip.H3K27me3 <- read.table('~/R/data/chipseq/hela/Hela-ChIP_H3K27me3_ENCFF000BBS_coverage_10000.bed', header = F, as.is = T)
colnames(Chip.H3K27me3) <- c('chr', 'start', 'end', 'reads_number', 'coverage', 'length', 'ratio')

Chip.H3K36me3 <- read.table('~/R/data/chipseq/hela/Hela-ChIP_H3K36me3_ENCFF000BCA_coverage_10000.bed', header = F, as.is = T)
colnames(Chip.H3K36me3) <- c('chr', 'start', 'end', 'reads_number', 'coverage', 'length', 'ratio')

Chip.H3K4me3 <- read.table('~/R/data/chipseq/hela/Hela-ChIP_H3K4me3_ENCFF000BCO_coverage_10000.bed', header = F, as.is = T)
colnames(Chip.H3K4me3) <- c('chr', 'start', 'end', 'reads_number', 'coverage', 'length', 'ratio')

Chip.H3K9me3 <- read.table('~/R/data/chipseq/hela/Hela-ChIP_H3K9me3_ENCFF000BBG_coverage_10000.bed', header = F, as.is = T)
colnames(Chip.H3K9me3) <- c('chr', 'start', 'end', 'reads_number', 'coverage', 'length', 'ratio')

Chip.H4K20me1 <- read.table('~/R/data/chipseq/hela/Hela-ChIP_H4K20me1_ENCFF000BDC_coverage_10000.bed', header = F, as.is = T)
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

plot.chip <- function(random.data, sample.data, range = c(50, 120), lwd  = 5, at.y = c(50, 75, 100),
                      col1 = '#FF9700', col2= '#2219b2', main = 'main', tck = -0.01){
  plot.new()
  plot.window(xlim = c(0, 110), ylim = range  )
  title(main = main)
  lines(1:101, apply(random.data, 2, score), col = col1, lwd = lwd * 0.5)
  lines(1:101, apply(sample.data, 2, score), col = col2, lwd = lwd * 0.67)
  axis(side = 1, lwd = 2,at = c(0, 51, 100), tck= tck)
  axis(side = 2, lwd = 2, at = at.y, tck= tck)
}

# NAIR plot 1 -------
#pdf('~/R/data/nucleulos/figures/NAIR_boundary_chip.pdf')
#png('~/R/data/nucleulos/figures/NAIR_boundary_chip.png')
layout(matrix(c(1,2,3,4, 5, 6), 2,3))
par(oma = c(4 , 4, 4 ,4))
par(mar = c(2, 0.5, 2, 0.5))

plot.chip(random.NAIR.boundary.Chip.H3K4me3, NAIR.boundary.Chip.H3K4me3, main = 'H3K4me3', range = c(50, 110), col2 = '#FB3F51', at.y = c(50, 75, 100))
plot.chip(random.NAIR.boundary.Chip.H3K9me3, NAIR.boundary.Chip.H3K9me3, range = c(100, 350), main = 'H3K9me3', col2 = '#D235D2', at.y = c(100, 200, 300))
plot.chip(random.NAIR.boundary.Chip.H4K20me1, NAIR.boundary.Chip.H4K20me1, main = 'H4K20me1', col2 = '#36D695', at.y = c(50, 80, 110))
plot.chip(random.NAIR.boundary.Chip.H3K27me3, NAIR.boundary.Chip.H3K27me3, range = c(100, 190), main = 'H3K27me3', col2 = '#33CCCC', at.y = c(100, 140, 180))
plot.chip(random.NAIR.boundary.Chip.H3K36me3, NAIR.boundary.Chip.H3K36me3, range = c(70, 150), main = 'H3K36me3', col2 = '#AC3BD4', at.y = c(80, 110, 140))
plot.chip(random.NAIR.boundary.Chip.CTCF, NAIR.boundary.Chip.CTCF, range = c(60, 120),main = 'CTCF', at.y = c(60, 90, 120))

#dev.off()

# NAIR plot 2 --------
plot.chip2 <- function(random.data, sample.data, range = c(50, 120), lwd  = 3, at.y = c(50, 75, 100),
                      col1 = '#FF9700', col2= '#2219b2', main = 'main', tck = -0.01){
  plot.new()
  plot.window(xlim = c(0, 110), ylim = range  )
  lines(1:101, apply(random.data, 2, score), col = col1, lwd = lwd * 0.34)
  lines(1:101, apply(sample.data, 2, score), col = col2, lwd = lwd * 0.67)
  axis(side = 1, lwd = 2,at = c(0, 51, 100), tck= tck, labels = F)
  axis(side = 2, lwd = 2, at = at.y, tck= tck, labels = F)
}

# NAIR plot
# pdf('~/R/data/nucleulos/figures/NAIR_boundary_chip2.pdf')
#png('~/R/data/nucleulos/figures/NAIR_boundary_chip.png')
layout(matrix(c(1,2,3,4, 5, 6), 2,3))
par(oma = c(4 , 4, 4 ,4))
par(mar = c(2, 0.5, 2, 0.5))

plot.chip2(random.NAIR.boundary.Chip.H3K4me3, NAIR.boundary.Chip.H3K4me3, main = 'H3K4me3', range = c(50, 110), col2 = '#FB3F51', at.y = c(50, 75, 100))
plot.chip2(random.NAIR.boundary.Chip.H3K9me3, NAIR.boundary.Chip.H3K9me3, range = c(100, 350), main = 'H3K9me3', col2 = '#D235D2', at.y = c(100, 200, 300))
plot.chip2(random.NAIR.boundary.Chip.H4K20me1, NAIR.boundary.Chip.H4K20me1, main = 'H4K20me1', col2 = '#36D695', at.y = c(50, 80, 110))
plot.chip2(random.NAIR.boundary.Chip.H3K27me3, NAIR.boundary.Chip.H3K27me3, range = c(100, 190), main = 'H3K27me3', col2 = '#33CCCC', at.y = c(100, 140, 180))
plot.chip2(random.NAIR.boundary.Chip.H3K36me3, NAIR.boundary.Chip.H3K36me3, range = c(70, 150), main = 'H3K36me3', col2 = '#AC3BD4', at.y = c(80, 110, 140))
plot.chip2(random.NAIR.boundary.Chip.CTCF, NAIR.boundary.Chip.CTCF, range = c(60, 120),main = 'CTCF', at.y = c(60, 90, 120))

#dev.off()


layout(matrix(1))
plot.new()
# TAD plot
plot.chip(random.NAIR.boundary.Chip.CTCF, TAD_boundary_Chip.CTCF)
plot.chip(random.NAIR.boundary.Chip.H3K4me3, TAD_boundary_Chip.H3K4me3)
plot.chip(random.NAIR.boundary.Chip.H3K9me3, TAD_boundary_Chip.H3K9me3, range = c(100, 200))
plot.chip(random.NAIR.boundary.Chip.H4K20me1, TAD_boundary_Chip.H4K20me1)
plot.chip(random.NAIR.boundary.Chip.H3K27me3, TAD_boundary_Chip.H3K27me3, range = c(100, 200))
plot.chip(random.NAIR.boundary.Chip.H3K36me3, TAD_boundary_Chip.H3K36me3, range = c(70, 150))



plot.window(xlim = c(0, 110), ylim =c(63, 116) )
title(main = main)
lines(1:101, apply(random_NAIR_boundary_Chip.CTCF, 2, score), col = '#FF9700', lwd = 5)
lines(1:101, apply(NAIR.boundary_Chip.CTCF, 2, score), col = '#2219b2', lwd = 5)

axis(side = 1, lwd = 2, at = c(0, 51, 100))
axis(side = 2, lwd = 2)


