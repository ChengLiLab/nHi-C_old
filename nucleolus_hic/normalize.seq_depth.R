
##############################################################################################
# calculate Hi-C experiment total interaction
##############################################################################################

# rm(list=ls())
# case_path = "~/menghw_HD/our_data/Hela-genome-hic-hg19/matrix/raw_matrix/MAPQ20/binsize_1000000/"
# binsize = 1e6
# matrix.path.format.cis = paste0(case_path,"chr_%d_%d_MAPQ20.txt")
# matrix.path.format.trans = paste0(case_path,"chr_%d_%d_%d_MAPQ20.txt")
# 
# 
# genome_matrix.sum = 0
# chrom_matrix.ratio = rep(0,23)
# 
# for(chrom_index.i in c(1:23)){
#   print(sprintf("Starting the Chr%d",chrom_index.i))
#   for(chrom_index.j in c(chrom_index.i:23)){
#     if(chrom_index.i == chrom_index.j){
#       matrix.path = sprintf(matrix.path.format.cis,chrom_index.i,binsize)
#       matrix.raw = read.table(matrix.path,header = F,sep = ",")
#       genome_matrix.sum = genome_matrix.sum + sum(as.matrix(matrix.raw))
#       # chrom_matrix.ratio[chrom_index.i] = sum(as.matrix(matrix.raw)) / genome_matrix.sum.hela_n_hic.raw
#       
#     }else{
#       matrix.path = sprintf(matrix.path.format.trans,chrom_index.i,chrom_index.j,binsize)
#       matrix.raw = read.table(matrix.path,header = F,sep = ",")
#       genome_matrix.sum = genome_matrix.sum + sum(as.matrix(matrix.raw)) * 2
#     }
#   }
# }


##############################################################################################
# normalize matrix and depart into cis matrix
##############################################################################################

#########################################
# record genome all heatmap interaction 
#########################################

# genome_matrix.sum.hela_g_hic.raw = 468290612
# genome_matrix.sum.hela_n_hic.raw = 499322076
# genome_matrix.sum.hela_n_hic_ActD.raw = 532089059
# 
# chrom_matrix.ratio.hela_n_hic_ActD.raw = c(0.048012904,0.033810879,0.025132041,0.017985049,0.039228707,
#                                            0.023333075,0.024004258,0.027140269,0.016936452,0.022072916,
#                                            0.020824903,0.021677341,0.010503418,0.008201753,0.011159707,
#                                            0.013617277,0.011318269,0.006808667,0.007247713,0.008663198,
#                                            0.003824798,0.006173085,0.017765567)
# 
# 
# chrom_matrix.ratio.hela_n_hic.raw = c(0.055666816,0.044325268,0.041199174,0.032095260,0.073684185,
#                                       0.038037912,0.039469194,0.048775534,0.028025302,0.034638012,
#                                       0.032060631,0.040551918,0.025068133,0.023886414,0.035000549,
#                                       0.026023896,0.016014804,0.020976717,0.012989139,0.016551147,
#                                       0.018815675,0.006142486,0.040188499)
# 
# chrom_matrix.ratio.hela_g_hic.raw = c(0.073356793,0.052524375,0.044827986,0.033402228,0.065420590,
#                                       0.043557465,0.039594712,0.043918505,0.028748308,0.034713076,
#                                       0.033996486,0.039560065,0.021658852,0.015304082,0.020546382,
#                                       0.019950708,0.018461691,0.013898214,0.010755591,0.013469016,
#                                       0.007758431,0.007725975,0.030270572)

# sequence_depth.hela_g_seq = 187807840
# sequence_depth.hela_n_seq = 126900190
# sequence_depth.hela_n_seq_ActD = 282845496
# 
# hic_count.hela_g_hic.trans = 55999521
# hic_count.hela_g_hic.cis = 247548509
# hic_count.hela_g_hic.total = hic_count.hela_g_hic.trans + hic_count.hela_g_hic.cis
# 
# hic_count.hela_n_hic.cis = 273219336
# hic_count.hela_n_hic.trans = 62986764
# hic_count.hela_n_hic.total = hic_count.hela_n_hic.trans + hic_count.hela_n_hic.cis
# 
# hic_count.hela_n_hic_ActD.cis = 137755064
# hic_count.hela_n_hic_ActD.trans = 153207862
# hic_count.hela_n_hic_ActD.total = hic_count.hela_n_hic_ActD.trans + hic_count.hela_n_hic_ActD.cis
# 
# save.image(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")



#########################################
# 矫正matrix的思路
# 1.计算频率矩阵
# 2.乘上系数
#########################################






##############################################################################################
# normalize sequencing depth 
##############################################################################################
# 












