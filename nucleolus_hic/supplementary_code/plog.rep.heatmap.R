########################################################################
# plot hic repeat heatmap
########################################################################
load(file = "~/menghw_HD/R_code/my_function/MyPlot_v08.RData")
hela_g_hic_rep1 = read.table(file = "~/menghw_HD/our_data/Hela-genome-hic-hg19-rep1/matrix/raw_matrix/chr_1_1000000_MAPQ20.txt",header = F,sep = ",")
hela_g_hic_rep2 = read.table(file = "~/menghw_HD/our_data/Hela-genome-hic-hg19-rep2/matrix/raw_matrix/chr_1_1000000_MAPQ20.txt",header = F,sep = ",")

par(mar=c(1,1,1,1))
plot.matrix(hela_g_hic_rep1,bound.max = 0.95,n_block_color = "darkgray")
plot.matrix(hela_g_hic_rep2,bound.max = 0.95,n_block_color = "darkgray")


hela_n_hic_rep1 = read.table(file = "~/menghw_HD/our_data/Hela-nucleolus-hic-hg19-rep1/matrix/raw_matrix/chr_1_1000000_MAPQ20.txt",header = F,sep = ",")
hela_n_hic_rep2 = read.table(file = "~/menghw_HD/our_data/Hela-nucleolus-hic-hg19-rep2/matrix/raw_matrix/chr_1_1000000_MAPQ20.txt",header = F,sep = ",")

par(mar=c(1,1,1,1))
plot.matrix(hela_n_hic_rep1,bound.max = 0.95,n_block_color = "darkgray")
plot.matrix(hela_n_hic_rep2,bound.max = 0.95,n_block_color = "darkgray")



