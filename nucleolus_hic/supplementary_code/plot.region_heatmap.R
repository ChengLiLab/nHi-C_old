##################################################################################
# plot region heatmap in fig5
##################################################################################
rm(list=ls())

################################################################
# plot matrix in one chromosome
################################################################
rm(list=ls())
chrom_index = 19
binsize = 50e3

load(file=sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_all_matrix_chr_%d_%d.RData",chrom_index,binsize))
load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")
hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
colnames(hela_NAIR_table) = c("chrom_name","start","end","peak_index","info")

chrom_name = chrom.name(chrom_index)
# conver as.matrix
segment.row = c(30e6,59e6)

# segment.row = c(0,chrom.length(chrom_index))
segment.col = segment.row
chrom_start = segment.col[1]
chrom_end = segment.col[2]

g_hic_matrix.part = matrix.part.corner(g_hic_matrix,chrom_index,segment.row,segment.col,binsize)
n_hic_matrix.part = matrix.part.corner(n_hic_matrix,chrom_index,segment.row,segment.col,binsize)
g_hic_ActD_matrix.part = matrix.part.corner(g_hic_ActD_matrix,chrom_index,segment.row,segment.col,binsize)
n_hic_ActD_matrix.part = matrix.part.corner(n_hic_ActD_matrix,chrom_index,segment.row,segment.col,binsize)


# make matrix with the same sequencing depth
g_hic_matrix.fix = g_hic_matrix.part / hic_count.hela_g_hic.total * 300e6
n_hic_matrix.fix = n_hic_matrix.part / hic_count.hela_n_hic.total * 300e6
g_hic_ActD_matrix.fix = g_hic_ActD_matrix.part / hic_count.hela_g_hic_ActD.total * 300e6
n_hic_ActD_matrix.fix = n_hic_ActD_matrix.part / hic_count.hela_n_hic_ActD.total * 300e6

quantile(g_hic_matrix.fix,prob=seq(0.9,1,0.01))
quantile(n_hic_matrix.fix,prob=seq(0.9,1,0.01))
quantile(g_hic_ActD_matrix.fix,prob=seq(0.9,1,0.01))
quantile(n_hic_ActD_matrix.fix,prob=seq(0.9,1,0.01))

par(mar=c(1,1,1,1))
plot.matrix(g_hic_matrix.fix,col.min = "red",col.max = "red",bound.max = 0.94,col.boundary = 0,n_block_color = "gray")
plot.matrix(n_hic_matrix.fix,col.min = "red",col.max = "red",bound.max = 0.95,col.boundary = 0,n_block_color = "gray")
plot.matrix(g_hic_ActD_matrix.fix,col.min = "red",col.max = "red",bound.max = 0.93,col.boundary = 0,n_block_color = "gray")
plot.matrix(n_hic_ActD_matrix.fix,col.min = "red",col.max = "red",bound.max = 0.96,col.boundary = 0,n_block_color = "gray")




