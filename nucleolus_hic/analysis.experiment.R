# analysis experiment enrichment effect

# example of chrom_index=16 enrichment region
# chrom_index = 16
# hela_g_hic_1mb = read.table("~/menghw_HD/our_data/Hela-genome-hic-hg19/matrix/raw_matrix/MAPQ20/binsize_1000000/chr_16_1000000_MAPQ20.txt",header = F,sep = ",")
# hela_g_hic_100kb = read.table("~/menghw_HD/our_data/Hela-genome-hic-hg19/matrix/raw_matrix/MAPQ20/binsize_100000/chr_16_100000_MAPQ20.txt",header = F,sep = ",")
# hela_g_hic_10kb = read.table("~/menghw_HD/our_data/Hela-genome-hic-hg19/matrix/raw_matrix/MAPQ20/binsize_20000/chr_16_20000_MAPQ20.txt",header = F,sep = ",")
# 
# hela_n_hic_1mb = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-hg19/matrix/raw_matrix/MAPQ20/binsize_1000000/chr_16_1000000_MAPQ20.txt",header = F,sep = ",")
# hela_n_hic_100kb = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-hg19/matrix/raw_matrix/MAPQ20/binsize_100000/chr_16_100000_MAPQ20.txt",header = F,sep = ",")
# hela_n_hic_10kb = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-hg19/matrix/raw_matrix/MAPQ20/binsize_20000/chr_16_20000_MAPQ20.txt",header = F,sep = ",")
# 
# hela_n_hic_ActD_1mb = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_1000000/chr_16_1000000_MAPQ20.txt",header = F,sep = ",")
# hela_n_hic_ActD_100kb = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_100000/chr_16_100000_MAPQ20.txt",header = F,sep = ",")
# hela_n_hic_ActD_10kb = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_20000/chr_16_20000_MAPQ20.txt",header = F,sep = ",")
# 
# save.image("~/menghw_HD/R_code/nucleolus_hic/r_image/chr_16_matrix.RData")

hg19_g_band_table = read.table("~/menghw_HD/reference/hg19_G-band.txt",header = F,sep = "\t")

load("~/menghw_HD/R_code/nucleolus_hic/r_image/chr_16_matrix.RData")
load("~/menghw_HD/R_code/my_function/MyPlot_v01.RData")

chrom_index = 16
chrom_start = 0
chrom_end = chrom.length(chrom_index)

g_hic_HiTCexp.obj = make_HTCexp_obj(hela_g_hic_10kb,chr_index_1 = chrom_index,chr_index_2 = chrom_index,resolution = 20e3)
g_hic_matrix.norm.obj = normICE(g_hic_HiTCexp.obj, 100)
g_hic_matrix.norm = as.matrix(intdata(g_hic_matrix.norm.obj))

png("~/menghw_HD/R_code/nucleolus_hic/r_figure/Hela-g-hic_chr16_100Kb")
plot.matrix(g_hic_matrix.norm,min_bound = 0,max_bound = 0.95)
print("done!")

# part enrichment region
chrom_start = 5e6
chrom_end = 5.5e6
binsize = 20e3
hela_g_hic_part = hela_g_hic_10kb[1150:1400,1150:1400]
hela_n_hic_part = hela_n_hic_10kb[1150:1400,1150:1400]

plot.matrix(hela_g_hic_part,min_bound = 0,max_bound = 0.98)
plot.matrix(hela_n_hic_part,min_bound = 0,max_bound = 0.90)

plot.peak(hg19_g_band_table,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom.length(chrom_index),color_pattern = T,track_color = "black",track_pattern = T,xaxt = F,yaxt = F,rect_border = T,rect_border.lwd = 2,xylwd = 3,xcex = 3)
par(new=T)
rect(xleft = 5e6,xright = 5.5e6,ybottom = -5,ytop = 5 ,border = T,col = "blue")
dev.off()
hg19_g_band_table[hg19_g_band_table$V1=="chr16",]


30e6 / 20e3

