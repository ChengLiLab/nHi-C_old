##################################################################################################
# rDNA matrix
##################################################################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v07.RData")
hela_g_hic_ActD_matrix.raw = read.table("~/menghw_HD/our_data/rDNA_matrix_hg19/Hela-g-hic-ActD_chr_25_500.txt",header = F,sep = ",")
hela_g_hic_matrix.raw = read.table("~/menghw_HD/our_data/rDNA_matrix_hg19/Hela-g-hic_chr_25_500.txt",header = F,sep = ",")

dev.off()
plot.matrix(hela_g_hic_matrix.raw,bound.max = 0.9)

plot.matrix(hela_g_hic_ActD_matrix.raw,bound.max = 0.9)


