# # # ###############################################################################################
# # # NAD region
# # # filter rule
# # # value > q(0.85) then merge within a 10Kb window size
# # ## hela-n-seq
# # ### broad peak 的效果不如 narrow peak的效果
# # hela_n_seq_broad.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-seq/hela-n-seq-broad-v1_peaks.broadPeak",header = F,sep = "\t")
# # hela_n_seq_narrow.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-seq/hela-n-seq-narrow-v1_peaks.narrowPeak",header = F,sep = "\t")
# # 
# # sum(hela_n_seq_broad.raw$V3 - hela_n_seq_broad.raw$V2)
# # sum(hela_n_seq_narrow.raw$V3 - hela_n_seq_narrow.raw$V2)
# # 
# # # quantile(hela_n_seq_broad.raw$V7)
# # # hela_n_seq_broad.fix = hela_n_seq_broad.raw[hela_n_seq_broad.raw$V7>=quantile(hela_n_seq_broad.raw$V7,prob=0.95),]
# # # sum(hela_n_seq_broad.fix$V3 - hela_n_seq_broad.fix$V2)
# # 
# # 
# # hela_n_seq_narrow.fix = hela_n_seq_narrow.raw[hela_n_seq_narrow.raw$V7>=quantile(hela_n_seq_narrow.raw$V7,prob=0.95),]
# # sum(hela_n_seq_narrow.fix$V3 - hela_n_seq_narrow.fix$V2)
# # nrow(hela_n_seq_narrow.fix)
# # write.table(hela_n_seq_narrow.fix,file = "~/menghw_HD/data_table/raw_table/Hela-n-seq_narrow_fix.table.v3",col.names = F,row.names = F,quote = F,sep = "\t")
# # 
# # 
# # # NAD_table = read.table(file = "~/menghw_HD/data_table/raw_table/Hela-n-seq_narrow_fix_w1e4.table")
# # # NAD_table = read.table(file = "~/menghw_HD/data_table/raw_table/Hela-n-seq_narrow_fix_w1e5.table")
# # NAD_table = read.table(file = "~/menghw_HD/data_table/raw_table/Hela-n-seq_narrow_fix_w1e5.table.v1")
# # hist((NAD_table$V3- NAD_table$V2))
# # quantile((NAD_table$V3- NAD_table$V2))
# # sum(NAD_table$V3- NAD_table$V2)
# # 
# # NAD_table.fix = NAD_table[(NAD_table$V3- NAD_table$V2)>= mean((NAD_table$V3- NAD_table$V2)), ]
# # 
# # quantile(NAD_table.fix$V3 - NAD_table.fix$V2)
# # sum(NAD_table.fix$V3 - NAD_table.fix$V2)
# # 
# # 
# # # 经过比较得到目前的NAD比较合理，合适
# # write.table(NAD_table.fix,file="~/menghw_HD/data_table/fix_table/Hela_NAD_region.bed",col.names = F,row.names = F,sep = "\t",quote = F)
# # 
# # 
# # # NAIR
# # rm(list=ls())
# # load(file="~/menghw_HD/R_code/my_function/MyPlot_v03.RData")
# # 
# # NAT_table.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-hic/seq-hic-control/hela-n-hic-v4_peaks.broadPeak")
# # NAT_table.raw.fix = NAT_table.raw[NAT_table.raw$V7>=quantile(NAT_table.raw$V7,prob=0.75),]
# # sum(NAT_table.raw.fix$V3 - NAT_table.raw.fix$V2)
# # 
# # write.table(NAT_table.raw.fix,file = "~/menghw_HD/data_table/raw_table/Hela-n-hic_seq_control_fix.table.v2",col.names = F,row.names = F,sep = "\t",quote = F)
# # 
# # 
# # NAT_table = read.table("~/menghw_HD/data_table/raw_table/Hela-n-hic_seq_control_fix_w1e5.bed.v2")
# # sum(NAT_table$V3 - NAT_table$V2)
# # quantile((NAT_table$V3 - NAT_table$V2))
# # 
# # NAT_table.fix = NAT_table[(NAT_table$V3-NAT_table$V2)>=median(NAT_table$V3-NAT_table$V2),]
# # sum(NAT_table.fix$V3 - NAT_table.fix$V2)
# # quantile((NAT_table.fix$V3 - NAT_table.fix$V2))
# # 
# # chrom_index = 1
# # chrom_start = 0
# # chrom_end = chrom.length(chrom_index)
# # 
# # png(file="~/menghw_HD/test5.png",width = 6000,height = 3000)
# # 
# # fig.facet <- layout(matrix(c(1:46),nrow = 46,byrow = FALSE),width = rep(10,46),heights = rep(1,46))
# # layout.show(fig.facet)
# # 
# # for(chrom_index in c(1:23)){
# #   par(mar=c(0.5,1,0,1))
# #   plot.peak(NAT_table.fix,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,
# #             yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#EE5C42",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
# #   par(mar=c(0.5,1,0,1))
# #   plot.chromosome(chrom_index,chrom_start,chrom_end)  
# # }
# # 
# # dev.off()
# # 
# # 
# # 
# # NAD_table = read.table("~/menghw_HD/data_table/fix_table/Hela_NAD_region.bed",header = F,sep = "\t")
# # 
# # write.table(NAT_table.fix,"~/menghw_HD/data_table/fix_table/Hela_NAT_region.bed",col.names = F,row.names = F,quote = F,sep = "\t")
# # 
# # NAD_NAT_overlap = overlap(NAD_table,NAT_table.fix,level_vetor = c(paste("chr",c(1:22),sep = ""),"chrX"))
# # sum(NAD_NAT_overlap$end - NAD_NAT_overlap$start)
# # 
# # 
# # quantile(hela_n_seq_narrow.raw$V5,prob=seq(0,1,0.1))
# # quantile_filter_value = quantile(hela_n_seq_narrow.raw$V5,prob=0.95)
# # quantile_filter_length = quantile(hela_n_seq_narrow.raw$V3 - hela_n_seq_narrow.raw$V2,prob=0.95)
# # hela_n_seq_narrow.fix = hela_n_seq_narrow.raw[hela_n_seq_narrow.raw$V5>=quantile_filter & (hela_n_seq_narrow.raw$V3 - hela_n_seq_narrow.raw$V2) >= quantile_filter_length,]
# # sum(hela_n_seq_narrow.fix$V3 - hela_n_seq_narrow.fix$V2)
# # quantile(hela_n_seq_narrow.fix$V3 - hela_n_seq_narrow.fix$V2)
# # 
# # write.table(hela_n_seq_narrow.fix,file="~/menghw_HD/data_table/raw_table/Hela-n-seq_narrow_fix.table",col.names = F,row.names = F,quote = F,sep = "\t")
# # 
# # hela_n_seq_narrow.fix.merge = read.table("~/menghw_HD/data_table/raw_table/Hela-n-seq_narrow_fix_w1e4.table",header = F,sep = "\t")
# # hela_n_seq_narrow.fix.merge = read.table("~/menghw_HD/data_table/raw_table/Hela-n-seq_narrow_fix_w1e5.table",header = F,sep = "\t")
# # 
# # sum(hela_n_seq_narrow.fix.merge$V3 - hela_n_seq_narrow.fix.merge$V2)
# # quantile((hela_n_seq_narrow.fix.merge$V3 - hela_n_seq_narrow.fix.merge$V2))
# # 
# # colnames(hela_n_seq_narrow.fix.merge) = c("chrom_name","start","end","name","value")
# # write.table(hela_n_seq_narrow.fix.merge,file="~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_w1e4.table",col.names = T,row.names = F,quote = F,sep = "\t")
# # 
# # ## hela-n-seq-ActD
# # hela_n_seq_ActD_narrow.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-seq-ActD/hela-n-seq-ActD-narrow-v1_peaks.narrowPeak",header = F,sep = "\t")
# # quantile_filter_value = quantile(hela_n_seq_ActD_narrow.raw$V5,prob=0.95)
# # quantile_filter_length = quantile(hela_n_seq_ActD_narrow.raw$V3 - hela_n_seq_ActD_narrow.raw$V2,prob=0.95)
# # 
# # hela_n_seq_ActD_narrow.fix = hela_n_seq_ActD_narrow.raw[hela_n_seq_ActD_narrow.raw$V5>=quantile_filter_value & (hela_n_seq_ActD_narrow.raw$V3 - hela_n_seq_ActD_narrow.raw$V2)>= quantile_filter_length,]
# # sum(hela_n_seq_ActD_narrow.fix$V3 - hela_n_seq_ActD_narrow.fix$V2)
# # write.table(hela_n_seq_ActD_narrow.fix,file="~/menghw_HD/data_table/raw_table/Hela-n-seq-ActD_narrow.fix",col.names = F,row.names = F,quote = F,sep = "\t")
# # 
# # hela_n_seq_ActD_narrow.fix.merge = read.table("~/menghw_HD/data_table/raw_table/Hela-n-seq-ActD_narrow_fix_w1e4.table",header = F,sep = "\t")
# # sum(hela_n_seq_ActD_narrow.fix.merge$V3 - hela_n_seq_ActD_narrow.fix.merge$V2)
# # 
# # colnames(hela_n_seq_ActD_narrow.fix.merge) = c("chrom_name","start","end","name","value")
# # write.table(hela_n_seq_narrow.fix.merge,file="~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_w1e4.table",col.names = T,row.names = F,quote = F,sep = "\t")
# # 
# # 
# # # hela-g-hic
# # ## 作为对照 与 hela-n-hic 的筛选条件相同
# # hela_g_hic_broad.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-g-hic/hela-g-hic-cis-broad-v1_peaks.broadPeak",header = F,sep = "\t")
# # 
# # sum(hela_g_hic_broad.raw$V3 - hela_g_hic_broad.raw$V2)
# # quantile(hela_g_hic_broad.raw$V5,prob=seq(0,1,0.1))
# # quantile((hela_g_hic_broad.raw$V3 - hela_g_hic_broad.raw$V2),prob=seq(0,1,0.1))
# # 
# # quantile_filter_value = quantile(hela_g_hic_broad.raw$V5,prob=0.85)
# # quantile_filter_length = quantile((hela_g_hic_broad.raw$V3 - hela_g_hic_broad.raw$V2),prob=0.5)
# # 
# # quantile_filter_value = 75
# # quantile_filter_length = 1000
# # 
# # hela_g_hic_broad.fix = hela_g_hic_broad.raw[hela_g_hic_broad.raw$V5>=quantile_filter_value &(hela_g_hic_broad.raw$V3 - hela_g_hic_broad.raw$V2)>= quantile_filter_length, ]
# # sum(hela_g_hic_broad.fix$V3 - hela_g_hic_broad.fix$V2)
# # quantile((hela_g_hic_broad.fix$V3 - hela_g_hic_broad.fix$V2))
# # 
# # write.table(hela_g_hic_broad.fix,file="~/menghw_HD/data_table/fix_table/Hela-g-hic_broad_fix.table",col.names = F,row.names = F,quote = F,sep = "\t")
# # 
# # hela_g_hic_broad.fix.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-g-hic_broad_fix_w1e5.table",header = F,sep = "\t")
# # sum(hela_g_hic_broad.fix.merge$V3 - hela_g_hic_broad.fix.merge$V2)
# # 
# # 
# # # hela-n-hic
# # hela_n_hic_broad.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-hic/seq-hic-control/hela-n-hic-cis-broad-v1_peaks.broadPeak",header = F,sep = "\t")
# # hela_n_hic_broad.length_vector = hela_n_hic_broad.raw$V3 - hela_n_hic_broad.raw$V2
# # sum(hela_n_hic_broad.length_vector)
# # quantile(hela_n_hic_broad.length_vector,prob=seq(0,1,0.1))
# # quantile(hela_n_hic_broad.raw$V5,prob=seq(0,1,0.1))
# # 
# # quantile_filter_value = 75
# # quantile_filter_length = 1000
# # 
# # hela_n_hic_broad.fix = hela_n_hic_broad.raw[hela_n_hic_broad.raw$V5 >= quantile_filter_value & hela_n_hic_broad.length_vector >= quantile_filter_length,]
# # sum(hela_n_hic_broad.fix$V3 - hela_n_hic_broad.fix$V2)
# # write.table(hela_n_hic_broad.fix,file="~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix.table",col.names = F,row.names = F,quote = F,sep = "\t")
# # 
# # hela_n_hic_broad.fix.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix_w1e5.table",header = F,sep = "\t")
# # sum(hela_n_hic_broad.fix.merge$V3 - hela_n_hic_broad.fix.merge$V2)
# # 
# # 
# # # hela-n-hic-ActD
# # hela_n_hic_ActD_broad.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-hic-ActD/Hela-n-hic-ActD-broad-v1_peaks.broadPeak",header = F,sep = "\t")
# # 
# # hela_n_hic_ActD_broad.length_vector = hela_n_hic_ActD_broad.raw$V3 - hela_n_hic_ActD_broad.raw$V2
# # sum(hela_n_hic_ActD_broad.length_vector)
# # 
# # quantile(hela_n_hic_ActD_broad.length_vector,prob=seq(0,1,0.1))
# # quantile(hela_n_hic_ActD_broad.raw$V5,prob=seq(0,1,0.1))
# # 
# # quantile_filter_value = 75
# # quantile_filter_length = 1000
# # 
# # hela_n_hic_ActD_broad.fix = hela_n_hic_ActD_broad.raw[hela_n_hic_ActD_broad.raw$V5 >= quantile_filter_value & hela_n_hic_ActD_broad.length_vector >= quantile_filter_length,]
# # sum(hela_n_hic_ActD_broad.fix$V3 - hela_n_hic_ActD_broad.fix$V2)
# # write.table(hela_n_hic_ActD_broad.fix,file="~/menghw_HD/data_table/fix_table/Hela-n-hic_ActD_broad_fix.table",col.names = F,row.names = F,quote = F,sep = "\t")
# # 
# # hela_n_hic_ActD_broad.fix.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_ActD_broad_fix_w1e5.table",header = F,sep = "\t")
# # sum(hela_n_hic_ActD_broad.fix.merge$V3 - hela_n_hic_ActD_broad.fix.merge$V2)
# # 
# # 
# # 
# # ############################################################################
# # # use fold change to identify the NAD and NAT
# # ############################################################################
# # 
# # rm(list=ls())
# # hela_all_table = read.table(file="~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep = "\t")
# # filter_vector = hela_all_table$not_n_rate >= 0.85
# # hela_all_table.fix = hela_all_table[filter_vector,]
# # binsize_kb = 100e3 / 1000
# # 
# # 
# # hela_g_seq_total = sum(hela_all_table.fix$g_count) / 1e6
# # g_seq_RPKM = hela_all_table.fix$g_count / hela_g_seq_total / binsize_kb * 2
# # 
# # hela_n_seq_total = sum(hela_all_table.fix$n_count) / 1e6
# # n_seq_RPKM = hela_all_table.fix$n_count / hela_n_seq_total / binsize_kb * 2
# # 
# # hela_g_hic_total = sum(hela_all_table.fix$g_hic) / 1e6
# # g_hic_RPKM = hela_all_table.fix$g_hic / hela_g_hic_total / binsize_kb * 2
# # 
# # hela_n_hic_total = sum(hela_all_table.fix$n_hic) / 1e6
# # n_hic_RPKM = hela_all_table.fix$n_hic / hela_n_hic_total / binsize_kb * 2
# # 
# # 
# # 
# # NAD_table = hela_all_table.fix[(n_seq_RPKM / g_seq_RPKM) >= quantile(n_seq_RPKM / g_seq_RPKM,prob=0.95,na.rm=T),]
# # nrow(NAD_table) 
# # NAD_table = NAD_table[!(is.na(NAD_table$region_min) | is.na(NAD_table$region_max)),]
# # 
# # NAT_table = hela_all_table.fix[(n_hic_RPKM / g_seq_RPKM) >=quantile(n_seq_RPKM / g_seq_RPKM,prob=0.95,na.rm=T), ]
# # NAT_table = NAT_table[!(is.na(NAT_table$region_min) | is.na(NAT_table$region_max)),]
# # 
# # 
# # load(file="~/menghw_HD/R_code/my_function/MyPlot_v03.RData")
# # 
# # chrom_index = 1
# # chrom_start = 0
# # chrom_end = chrom.length(chrom_index)
# # 
# # png(file="~/menghw_HD/test6.png",width = 6000,height = 3000)
# # 
# # fig.facet <- layout(matrix(c(1:46),nrow = 46,byrow = FALSE),width = rep(10,46),heights = rep(1,46))
# # layout.show(fig.facet)
# # 
# # for(chrom_index in c(1:23)){
# #   par(mar=c(0.5,1,0,1))
# #   plot.peak(NAT_table.fix,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,
# #             yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#EE5C42",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
# #   par(mar=c(0.5,1,0,1))
# #   plot.chromosome(chrom_index,chrom_start,chrom_end)  
# # }
# # 
# # dev.off()
# # 
# # write.table(NAD_table,file="~/menghw_HD/data_table/raw_table/Hela-n-seq_NAD_test.bed",col.names = F,row.names = F,quote = F,sep = "\t")
# # write.table(NAT_table,file="~/menghw_HD/data_table/raw_table/Hela-n-seq_NAT_test.bed",col.names = F,row.names = F,quote = F,sep = "\t")
# # 
# # NAD_table.fix = read.table("~/menghw_HD/data_table/raw_table/Hela-n-seq_NAD_w2e5.bed",header = F,sep = "\t")
# # NAT_table.fix = read.table("~/menghw_HD/data_table/raw_table/Hela-n-seq_NAT_w2e5.bed",header = F,sep = "\t")
# # 
# # NAD_NAT_overlap = overlap(NAD_table.fix,NAT_table.fix,level_vetor = c(paste("chr",c(1:22),sep = ""),"chrX"))
# # 
# # sum(NAD_table.fix$V3 - NAD_table.fix$V2)
# # quantile((NAD_table.fix$V3 - NAD_table.fix$V2))
# # 
# # sum(NAT_table.fix$V3 - NAT_table.fix$V2)
# # quantile(NAT_table.fix$V3 - NAT_table.fix$V2)
# # 
# # 
# # 
# # NAD_table.fix.fix = NAD_table.fix
# # NAD_table.fix.fix$V2 = NAD_table.fix$V2 + ceiling(rnorm(nrow(NAD_table.fix),mean = 0,sd=10e3))
# # NAD_table.fix.fix$V3 = NAD_table.fix$V3 + ceiling(rnorm(nrow(NAD_table.fix),mean = 0,sd=10e3))
# # quantile(NAD_table.fix.fix$V3 - NAD_table.fix.fix$V2,prob=seq(0,1,0.05))
# # NAD_table.fix.fix = NAD_table.fix.fix[(NAD_table.fix.fix$V3 - NAD_table.fix.fix$V2)>=quantile(NAD_table.fix.fix$V3 - NAD_table.fix.fix$V2,prob=0.1), ]
# # sum((NAD_table.fix.fix$V3 - NAD_table.fix.fix$V2))
# # 
# # 
# # 
# # NAT_table.fix.fix = NAT_table.fix
# # NAT_table.fix.fix$V2 = NAT_table.fix$V2 + ceiling(rnorm(nrow(NAT_table.fix),mean = 0,sd=10e3))
# # NAT_table.fix.fix$V3 = NAT_table.fix$V3 + ceiling(rnorm(nrow(NAT_table.fix),mean = 0,sd=10e3))
# # quantile(NAT_table.fix.fix$V3 - NAT_table.fix.fix$V2,prob=seq(0,1,0.1))
# # NAT_table.fix.fix = NAT_table.fix.fix[(NAT_table.fix.fix$V3 - NAT_table.fix.fix$V2)>=quantile(NAT_table.fix.fix$V3 - NAT_table.fix.fix$V2,prob=0.1), ]
# # sum((NAT_table.fix.fix$V3 - NAT_table.fix.fix$V2))
# # 
# # 
# # 
# # write.table(NAD_table.fix.fix,file = "~/menghw_HD/data_table/Hela_NAD.bed",col.names = F,row.names = F,quote = F,sep = "\t")
# # write.table(NAT_table.fix.fix,file = "~/menghw_HD/data_table/Hela_NAT.bed",col.names = F,row.names = F,quote = F,sep = "\t")
# # 
# # i = 1
# # for(i in c(1:23)){
# #   chrom_name = chrom.name(i)
# #   chrom_length = chrom.length(i)
# #   
# #   NAD_table.fix.fix[NAD_table.fix.fix$V1==chrom_name,]$V2[NAD_table.fix.fix[NAD_table.fix.fix$V1==chrom_name,]$V2 < 0] = 0
# #   NAD_table.fix.fix[NAD_table.fix.fix$V1==chrom_name,]$V3[NAD_table.fix.fix[NAD_table.fix.fix$V1==chrom_name,]$V3 > chrom_length] = chrom_length
# #   
# #   NAT_table.fix.fix[NAT_table.fix.fix$V1==chrom_name,]$V2[NAT_table.fix.fix[NAT_table.fix.fix$V1==chrom_name,]$V2 < 0] = 0
# #   NAT_table.fix.fix[NAT_table.fix.fix$V1==chrom_name,]$V3[NAT_table.fix.fix[NAT_table.fix.fix$V1==chrom_name,]$V3 > chrom_length] = chrom_length
# # }
# # 
# # 
# # 
# # NAD_NAT_overlap = overlap(NAD_table.fix.fix,NAT_table.fix.fix,level_vetor = c(paste("chr",c(1:22),sep = ""),"chrX"))
# # sum(NAD_NAT_overlap$end - NAD_NAT_overlap$start) / sum(NAD_table.fix.fix$V3 - NAD_table.fix.fix$V2)
# 
# ############################################################################
# # use fold change to identify the NAD and NAT, fold change cutoff = 2
# ############################################################################
# rm(list=ls())
# #################################
# # sequencing data
# #################################
# hela_n_seq_total = 126900190 * 135
# hela_n_seq_ActD_total = 282845496 * 135
# hela_g_seq_total = 187807841 * 135
# 
# hela_g_hic_trans = 55999521 * 2 * 36
# hela_g_hic_cis = 247548509 * 2 * 36
# 
# hela_n_hic_trans = 62986764 * 2 * 36
# hela_n_hic_cis = 273219336 * 2 * 36
# 
# hela_n_hic_ActD_trans = 153207862 * 2 * 36
# hela_n_hic_ActD_cis = 137755064 * 2 * 36
# 
# #################################
# # call NAD
# #################################
# hela_all_table = read.table(file="~/menghw_HD/data_table/Hela_all_table_10000.txt",header = T,sep = "\t")
# filter_vector_NAD = hela_all_table$not_n_rate >= 0.85 & hela_all_table$g_count >0 & hela_all_table$g_cover_rate >=0.5
# hela_all_table.fix.NAD = hela_all_table[filter_vector_NAD,]
# 
# ######################
# # hela-n-seq NAD
# ######################
# # fix_value = hela_g_seq_total / hela_n_seq_total
# fix_value = 1.12
# hela_NAD_ratio = hela_all_table.fix.NAD$n_count / hela_all_table.fix.NAD$g_count * fix_value
# quantile(hela_NAD_ratio,seq(0.5,1,0.05))
# hela_NAD_table.raw = hela_all_table.fix.NAD[hela_NAD_ratio>2,]
# 
# write.table(hela_NAD_table.raw[,c(1,2,3)],file = "~/menghw_HD/data_table/Hela_NAD_fold2.bed",col.names = F,row.names = F,quote = F,sep = "\t")
# 
# hela_NAD.fix = read.table(file = "~/menghw_HD/data_table/Hela_NAD_fold2_w1e5.bed",header = F,sep = '\t')
# length_vector = as.vector(hela_NAD.fix$V3 - hela_NAD.fix$V2)
# hela_NAD.fix.fix = hela_NAD.fix[length_vector>120e3,]
# sum(hela_NAD.fix.fix$V3 - hela_NAD.fix.fix$V2)
# hist(log10(hela_NAD.fix.fix$V3 - hela_NAD.fix.fix$V2))
# quantile((hela_NAD.fix.fix$V3 - hela_NAD.fix.fix$V2))
# 
# 
# hela_NAD.fix.fix$V3 =  hela_NAD.fix.fix$V3 + rnorm(nrow(hela_NAD.fix.fix),mean = 1000,sd = 600)
# hela_NAD.fix.fix$V2 =  hela_NAD.fix.fix$V2 + rnorm(nrow(hela_NAD.fix.fix),mean = 1000,sd = 600)
# sum(hela_NAD.fix.fix$V3 - hela_NAD.fix.fix$V2)
# hela_NAD.fix.fix$V3 = round(hela_NAD.fix.fix$V3)
# hela_NAD.fix.fix$V2 = round(hela_NAD.fix.fix$V2)
# 
# write.table(hela_NAD.fix.fix,file = "~/menghw_HD/data_table/Hela_NAD_v2.bed",col.names = F,row.names = F,quote = F,sep = "\t")
# 
# #######################
# # hela-n-seq-ActD NAD
# #######################
# fix_value = hela_g_seq_total / hela_n_seq_ActD_total
# fix_value = 0.6
# hela_ActD_NAD_ratio = hela_all_table.fix.NAD$n_actd_count / hela_all_table.fix.NAD$g_count * fix_value
# quantile(hela_ActD_NAD_ratio,seq(0.5,1,0.05))
# hela_NAD_ActD_table.raw = hela_all_table.fix.NAD[hela_ActD_NAD_ratio>=2,]
# 
# 
# 
# #################################
# # call NAT
# #################################
# filter_vector_NAT = hela_all_table$not_n_rate >= 0.85 & hela_all_table$g_hic >0 & hela_all_table$g_cover_rate >=0.5
# hela_all_table.fix.NAT = hela_all_table[filter_vector_NAT,]
# 
# fix_value = ( hela_g_hic_cis) / (hela_n_hic_cis ) 
# # fix_value = 1
# 
# hela_NAT_ratio = hela_all_table.fix.NAT$n_hic / hela_all_table.fix.NAT$g_hic * fix_value
# 
# quantile(hela_NAT_ratio,seq(0.5,1,0.05))
# hela_NAT_table.raw = hela_all_table.fix.NAT[hela_NAT_ratio>=1.5,]
# 
# write.table(hela_NAT_table.raw[,c(1,2,3)],file = "~/menghw_HD/data_table/Hela_NAT_fold2.bed",col.names = F,row.names = F,quote = F,sep = "\t")
# 
# hela_NAT.fix = read.table(file = "~/menghw_HD/data_table/Hela_NAT_fold2_w1e5.bed",header = F,sep = '\t')
# 
# length_vector = as.vector(hela_NAT.fix$V3 - hela_NAT.fix$V2)
# quantile(length_vector)
# hela_NAT.fix.fix = hela_NAT.fix[length_vector>120e3,]
# 
# sum(hela_NAT.fix.fix$V3 - hela_NAT.fix.fix$V2)
# 
# hist(log10(hela_NAT.fix.fix$V3 - hela_NAT.fix.fix$V2))
# quantile((hela_NAT.fix.fix$V3 - hela_NAT.fix.fix$V2))
# hela_NAT.fix.fix[(hela_NAT.fix.fix$V3 - hela_NAT.fix.fix$V2)>5e6,]
# 
# 
# hela_NAT.fix.fix$V3 =  hela_NAT.fix.fix$V3 + rnorm(nrow(hela_NAT.fix.fix),mean = 1000,sd = 600)
# hela_NAT.fix.fix$V2 =  hela_NAT.fix.fix$V2 + rnorm(nrow(hela_NAT.fix.fix),mean = 1000,sd = 600)
# sum(hela_NAT.fix.fix$V3 - hela_NAT.fix.fix$V2)
# hela_NAT.fix.fix$V3 = round(hela_NAT.fix.fix$V3)
# hela_NAT.fix.fix$V2 = round(hela_NAT.fix.fix$V2)
# 
# write.table(hela_NAD.fix.fix,file = "~/menghw_HD/data_table/Hela_NAT_v2.bed",col.names = F,row.names = F,quote = F,sep = "\t")
# 
# 
# ##################################################################
# # use png test B compartment and enrichment region
# ##################################################################
# load(file="~/menghw_HD/R_code/my_function/MyPlot_v05.RData")
# hela_NAD_table = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep = "\t")
# hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
# hela_ab_value_table = read.table("~/menghw_HD/data_table/Hela_AB_compart_value.table",header = T,sep = "\t")
# 
# binsize = 100e3
# chrom_index = 1
# chrom_name = chrom.name(chrom_index)
# chrom_start = 0
# chrom_end = chrom.length(chrom_index)
# 
# image_path = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_track_chr_%d_%d.RData",chrom_index,binsize)
# load(image_path)
# png_path = sprintf("~/menghw_HD/R_image/Hela-hic-figure/20160830/Hela_matrix_chr_%d_%d_v1.png",chrom_index,binsize)
# 
# png(png_path,width = 2000,height = 3000)
# fig.facet <- layout(matrix(c(1:6),nrow = 6,byrow = FALSE),width = rep(10,6),heights = c(10,1,1,1,10,1))
# layout.show(fig.facet)
# 
# ## g_hic_matrix.ice.norm
# par(mar=c(0.5,1,1,1))
# plot.matrix(matrix.part(g_hic_matrix.norm,chrom_index,chrom_start,chrom_end,binsize),max_bound = 0.96)
# # plot.matrix(g_hic_matrix.norm,max_bound = 0.95)
# # box(lwd=3)
# 
# ## A,B compartment
# par(mar=c(0.5,1,0,1))
# ab_x = start(pc.norm$PC1)
# ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
# plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",yaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
# 
# ## NAD
# par(mar=c(0.5,1,0,1))
# plot.peak(hela_NAD_table,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#EE5C42",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
# 
# ## NAT
# par(mar=c(0.5,1,0,1))
# plot.peak(hela_NAT_table,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom_len,yaxt = F,xaxt = F,track_height = 2,axis_x_len=axis_x_len,track_color = "#EE5C42",track_pattern = T,color_pattern = F,peak.lwd = 0.5,xylwd = 3)
# 
# ##n_hic_matrix
# par(mar=c(0.1,1,0,1))
# plot.matrix(matrix.part(n_hic_matrix,chrom_index,chrom_start,chrom_end,binsize),max_bound = 0.96)
# 
# ## chrom 坐标轴
# par(mar=c(0.5,1,0.5,1))
# plot.chromosome(chrom_index,chrom_start,chrom_end,border.lwd = 2)
# 
# dev.off()
# print(sprintf("The chr%d ice.norm figure is done!",chrom_index))




#############################################################################################################
# use fold change to identify the NAD and NAT, fold change cutoff = 2
#############################################################################################################
rm(list=ls())

##################################################################
# sequencing data
##################################################################
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")

##################################################################
# call NAD
##################################################################
# use 10 Kb hela_all_table
hela_all_table = read.table("~/menghw_HD/data_table/Hela_all_table_10000.txt",header = T,sep = "\t")

##########################################
# hela-n-seq NAD & hela-n-seq-ActD NAD
##########################################
hela_all_table.raw = hela_all_table
hela_all_table = hela_all_table[hela_all_table$not_n_rate>=0.5 & hela_all_table$g_cover_rate >= 0.5,]


NAD_ratio = (hela_all_table$n_count ) / (hela_all_table$g_count)
NAD_ActD_ratio = (hela_all_table$n_actd_count / sequence_depth.hela_n_seq_ActD) / (hela_all_table$g_count / sequence_depth.hela_g_seq)

NAD_ratio[hela_all_table$g_count == 0 ] = 0
NAD_ActD_ratio[hela_all_table$g_count == 0 ] = 0

length(NAD_ratio[NAD_ratio>2]) * 10e3 / 1e6
length(NAD_ActD_ratio[NAD_ActD_ratio>2]) * 10e3 / 1e6

hela_NAD_table = cbind(hela_all_table[NAD_ratio>2.5,c(1,2,3,4,6)],NAD_ratio[NAD_ratio>2.5])
colnames(hela_NAD_table) = c(colnames(hela_NAD_table)[1:5],"NAD_ratio")

hela_ActD_NAD_table = cbind(hela_all_table[NAD_ActD_ratio>2.5,c(1,2,3,4,6)],NAD_ActD_ratio[NAD_ActD_ratio>2.5])
colnames(hela_ActD_NAD_table) = c(colnames(hela_ActD_NAD_table)[1:5],"NAD_ActD_ratio")


write.table(hela_NAD_table,file="~/menghw_HD/data_table/Hela_NAD_v4_raw.table",col.names = F,row.names = F,sep = "\t",quote = F)
write.table(hela_ActD_NAD_table,file="~/menghw_HD/data_table/Hela_ActD_NAD_v4_raw.table",col.names = F,row.names = F,sep = "\t",quote = F)

##########################################
# NAD: filter the outlier too small
##########################################
rm(list=ls())
hela_NAD_table.raw = read.table("~/menghw_HD/data_table/enrichment_region/Hela_NAD_v4_w1e5.table",header = F,sep = "\t")
hela_ActD_NAD_table.raw = read.table("~/menghw_HD/data_table/enrichment_region/Hela_ActD_NAD_v4_w1e5.table",header = F,sep = "\t")

quantile(hela_NAD_table.raw$V3 - hela_NAD_table.raw$V2)
quantile(hela_ActD_NAD_table.raw$V3 - hela_ActD_NAD_table.raw$V2)

hela_NAD_table.filter = hela_NAD_table.raw[(hela_NAD_table.raw$V3 - hela_NAD_table.raw$V2) > 100e3,]
hela_NAD_table.filter$V2 = hela_NAD_table.filter$V2 - abs(ceiling(rnorm(nrow(hela_NAD_table.filter),5000,5000)))
hela_NAD_table.filter$V3 = hela_NAD_table.filter$V3 + abs(ceiling(rnorm(nrow(hela_NAD_table.filter),5000,5000)))

sum(hela_NAD_table.filter$V3 - hela_NAD_table.filter$V2)
quantile((hela_NAD_table.filter$V3 - hela_NAD_table.filter$V2))
quantile(hela_ActD_NAD_table.raw$V3 -hela_ActD_NAD_table.raw$V2)

write.table(hela_NAD_table.filter,file="~/menghw_HD/data_table/Hela_NAD.bed",col.names = F,row.names = F,sep = "\t",quote = F)

####################################################################
# call NAT
####################################################################
# use matrix with 20K binsize
rm(list=ls())

###################################
# hela-hic colsum data.frame 
###################################

load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")

binsize = 20e3
genome_region_table = NULL

for(chrom_index in c(1:23)){
  print(sprintf("Starting chr%d",chrom_index))
  image_path = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_%d_%d.RData",chrom_index,binsize)
  load(image_path)
  matrix_path = sprintf("~/menghw_HD/our_data/Hela-nucleolus-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_%d/chr_%d_%d_MAPQ20.txt",binsize,chrom_index,binsize)
  n_hic_ActD_matrix = read.table(matrix_path,header = F,sep = ",")
  
  g_hic_colsum = colSums(g_hic_matrix)
  n_hic_colsum = colSums(n_hic_matrix)
  n_hic_ActD_colsum = colSums(n_hic_ActD_matrix)
  
  chrom_name = chrom.name(chrom_index)
  chrom_len = chrom.length(chrom_index)
  region_start = seq(0,chrom_len %/% binsize * binsize,binsize)
  region_end = region_start + binsize - 1
  chrom_region_table = data.frame(chrom_name = chrom_name,
                                  region_start = region_start,
                                  region_end = region_end,
                                  g_hic_colsum = g_hic_colsum,
                                  n_hic_colsum = n_hic_colsum,
                                  n_hic_ActD_colsum = n_hic_ActD_colsum)
  genome_region_table = rbind(genome_region_table,chrom_region_table)
}

genome_region_table.fix = genome_region_table
genome_region_table.fix$region_start = as.integer(genome_region_table$region_start)
genome_region_table.fix$region_end = as.integer(genome_region_table$region_end)

write.table(genome_region_table.fix,file = "~/menghw_HD/data_table/Hela_hic_colsum.table",col.names = T,row.names = F,sep = "\t",quote = F)

###################################
# hela-NAT
###################################
rm(list=ls())
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")
hic_colsum_table = read.table("~/menghw_HD/data_table/Hela_hic_colsum.table",header = T,sep = "\t")

NAT_ratio = (hic_colsum_table$n_hic_colsum / hic_count.hela_n_hic.total) / ((hic_colsum_table$g_hic_colsum) / hic_count.hela_g_hic.total)
NAT_ActD_ratio = (hic_colsum_table$n_hic_ActD_colsum / hic_count.hela_n_hic_ActD.total) / ((hic_colsum_table$g_hic_colsum) / hic_count.hela_g_hic.total)

NAT_ratio[hic_colsum_table$g_hic_colsum==0] = 0 
NAT_ActD_ratio[hic_colsum_table$g_hic_colsum==0] = 0

length(NAT_ratio[NAT_ratio>1.414]) * 20e3 / 1e6
length(NAT_ActD_ratio[NAT_ActD_ratio>1.414]) * 20e3 / 1e6


hela_NAT_table = cbind(hic_colsum_table[NAT_ratio>1.414,],NAT_ratio[NAT_ratio>1.414])
colnames(hela_NAT_table) = c(colnames(hela_NAT_table)[1:6],"NAT_ratio")

hela_ActD_NAT_table = cbind(hic_colsum_table[NAT_ActD_ratio>1.414,],NAT_ActD_ratio[NAT_ActD_ratio>1.414])
colnames(hela_ActD_NAT_table) = c(colnames(hela_ActD_NAT_table)[1:6],"NAT_ActD_ratio")


write.table(hela_NAT_table,file="~/menghw_HD/data_table/Hela_NAT_v4_raw.table",col.names = F,row.names = F,sep = "\t",quote = F)
write.table(hela_ActD_NAT_table,file="~/menghw_HD/data_table/Hela_ActD_NAT_v4_raw.table",col.names = F,row.names = F,sep = "\t",quote = F)


##########################################
# NAT: filter the outlier too small
##########################################
rm(list=ls())
hela_NAT_table.raw = read.table("~/menghw_HD/data_table/enrichment_region/Hela_NAT_v4_w1e5.table",header = F,sep = "\t")
hela_ActD_NAT_table.raw = read.table("~/menghw_HD/data_table/enrichment_region/Hela_ActD_NAT_v4_w1e5.table",header = F,sep = "\t")

quantile(hela_NAT_table.raw$V3 - hela_NAT_table.raw$V2)
quantile(hela_ActD_NAT_table.raw$V3 - hela_ActD_NAT_table.raw$V2)

hela_NAT_table.filter = hela_NAT_table.raw[(hela_NAT_table.raw$V3 - hela_NAT_table.raw$V2) > 100e3,]
hela_NAT_table.filter$V2 = hela_NAT_table.filter$V2 - abs(ceiling(rnorm(nrow(hela_NAT_table.filter),5000,5000)))
hela_NAT_table.filter$V3 = hela_NAT_table.filter$V3 + abs(ceiling(rnorm(nrow(hela_NAT_table.filter),5000,5000)))

sum(hela_NAT_table.filter$V3 - hela_NAT_table.filter$V2)
quantile((hela_NAT_table.filter$V3 - hela_NAT_table.filter$V2))
write.table(hela_NAT_table.filter,file="~/menghw_HD/data_table/Hela_NAT.bed",col.names = F,row.names = F,sep = "\t",quote = F)


hela_ActD_NAT_table.raw = read.table("~/menghw_HD/data_table/enrichment_region/Hela_ActD_NAT_v4_w1e5.table",header = F,sep = "\t")
quantile(hela_ActD_NAT_table.raw$V3 - hela_ActD_NAT_table.raw$V2)
sum((hela_ActD_NAT_table.raw$V3 - hela_ActD_NAT_table.raw$V2))

hela_ActD_NAT_table.filter = hela_ActD_NAT_table.raw[(hela_ActD_NAT_table.raw$V3 - hela_ActD_NAT_table.raw$V2) > 100e3,]
hela_ActD_NAT_table.filter$V2 = hela_ActD_NAT_table.filter$V2 - abs(ceiling(rnorm(nrow(hela_ActD_NAT_table.filter),5000,5000)))
hela_ActD_NAT_table.filter$V3 = hela_ActD_NAT_table.filter$V3 + abs(ceiling(rnorm(nrow(hela_ActD_NAT_table.filter),5000,5000)))
quantile(hela_ActD_NAT_table.filter$V3 -hela_ActD_NAT_table.filter$V2)

sum((hela_ActD_NAT_table.filter$V3 -hela_ActD_NAT_table.filter$V2))




###################################
# U2OS NAIR 区域
###################################
u2os_g_seq_table = read.table("~/menghw_HD/our_data/U2OS-genome-sequence-hg19/tmp.data/BAM/U2OS-genome-seq-hg19_MAPQ20_coverage_w1e5.bed",header = F,sep = "\t")
colnames(u2os_g_seq_table) = c("chrom_name","start","end","count","coverage","binsize","coverage_rate")

u2os_n_hic_table = read.table("~/menghw_HD/our_data/U2OS-nucleolus-hic-hg19/tmp.data/filtered/U2OS-n-hic_coverage_w1e5.bed",header = F,sep = "\t")
colnames(u2os_n_hic_table) = c("chrom_name","start","end","count","coverage","binsize","coverage_rate")

u2os_ratio= (u2os_n_hic_table$count / sum(u2os_n_hic_table$count)) / (u2os_g_seq_table$count / sum(u2os_g_seq_table$count))
u2os_ratio_table = cbind(u2os_n_hic_table[,c(1,2,3)],u2os_ratio)
colnames(u2os_ratio_table) = c("chrom_name","start","end","u2os_ratio")
u2os_ratio_table$u2os_ratio[u2os_g_seq_table$count==0] = 0

quantile(u2os_ratio_table$u2os_ratio)
hist(u2os_ratio_table$u2os_ratio[u2os_ratio_table$u2os_ratio>0 & u2os_ratio_table$u2os_ratio<5],xlim=c(0,5),breaks = seq(0,5,0.005))


u2os_NAIR_v1 = u2os_ratio_table.fix[u2os_ratio_table$u2os_ratio>= 1.414 & u2os_ratio_table$u2os_ratio <= 10,]
u2os_NAIR_v2 = u2os_ratio_table.fix[u2os_ratio_table$u2os_ratio>= 2 & u2os_ratio_table$u2os_ratio <= 10,]

write.table(u2os_NAIR_v1,file = "~/menghw_HD/data_table/u2os_NAIR_data/U2OS_NAIR_v1.bed",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(u2os_NAIR_v2,file = "~/menghw_HD/data_table/u2os_NAIR_data/U2OS_NAIR_v2.bed",col.names = T,row.names = F,quote = F,sep = "\t")

u2os_NAIR_v1.merge = read.table("~/menghw_HD/data_table/u2os_NAIR_data/U2OS_NAIR_v1_w1e5.bed")
u2os_NAIR_v1.merge.fix = u2os_NAIR_v1.merge
u2os_NAIR_v1.merge.fix$V2 = u2os_NAIR_v1.merge$V2 - abs(round(rnorm(length(u2os_NAIR_v1.merge$V2),5000,5000)))
u2os_NAIR_v1.merge.fix$V3 = u2os_NAIR_v1.merge$V3 +  abs(round(rnorm(length(u2os_NAIR_v1.merge$V2),5000,5000)))
u2os_NAIR_v1.merge.fix$V2[u2os_NAIR_v1.merge.fix$V2<0] = 0
colnames(u2os_NAIR_v1.merge.fix) = c("chrom_name","start","end","peak-index","test")
write.table(u2os_NAIR_v1.merge.fix,file = "~/menghw_HD/data_table/u2os_NAIR_data/U2OS_NAIR_v1_fix.bed",col.names = T,row.names = F,quote = F,sep = "\t")

u2os_NAIR_v2.merge = read.table("~/menghw_HD/data_table/u2os_NAIR_data/U2OS_NAIR_v2_w1e5.bed")
u2os_NAIR_v2.merge.fix = u2os_NAIR_v2.merge
u2os_NAIR_v2.merge.fix$V2 = u2os_NAIR_v2.merge$V2 - abs(round(rnorm(length(u2os_NAIR_v2.merge$V2),5000,5000)))
u2os_NAIR_v2.merge.fix$V3 = u2os_NAIR_v2.merge$V3 +  abs(round(rnorm(length(u2os_NAIR_v2.merge$V2),5000,5000)))
u2os_NAIR_v2.merge.fix$V2[u2os_NAIR_v2.merge.fix$V2<0] = 0
colnames(u2os_NAIR_v2.merge.fix) = c("chrom_name","start","end","peak-index","test")
write.table(u2os_NAIR_v2.merge.fix,file = "~/menghw_HD/data_table/u2os_NAIR_data/U2OS_NAIR_v2_fix.bed",col.names = T,row.names = F,quote = F,sep = "\t")

u2os_NAIR_v1 = read.table("~/menghw_HD/data_table/u2os_NAIR_data/U2OS_NAIR_v1_fix.bed",header = T,sep = '\t')
sum(u2os_NAIR_v1$end - u2os_NAIR_v1$start)
quantile(u2os_NAIR_v1$end - u2os_NAIR_v1$start)



