# Howard MENG on 2016-July-1st
# make enrichment region statistics
# We deal peak calling with MACS2 software combine with --broad parameter

################################################################################################
# constent parameters
# hg19 长度
HG19_LEN.total = 3095720411
HG19_LEN.n = 234350292
HG19_LEN.not_n = 2861370119

################################################################################################
# enrichment table 
# make LAD region 
# cLAD_table_raw = read.table("~/menghw_HD/data_table/public_table/GSE22428_LAD_overlap_w1e4.txt",header = F,sep = "\t")
# quantile((cLAD_table_raw$V3 - cLAD_table_raw$V2),prob=seq(0,1,0.1))
# quantile_filter = quantile((cLAD_table_raw$V3 - cLAD_table_raw$V2),prob=0.85)
# cLAD_table_fix = cLAD_table_raw[(cLAD_table_raw$V3 - cLAD_table_raw$V2)>=quantile_filter,]
# sum(cLAD_table_fix$V3 - cLAD_table_fix$V2)
# quantile((cLAD_table_fix$V3 - cLAD_table_fix$V2),prob=seq(0,1,0.1))
# cLAD_table_raw[(cLAD_table_raw$V3 - cLAD_table_raw$V2) > 1e6,]
# 
# colnames(cLAD_table_fix) = c("chrom_name","start","end","name","value")
# write.table(cLAD_table_fix,file="~/menghw_HD/data_table/fix_table/cLAD_human_fix.table",col.names = T,row.names = F,quote = F,sep = "\t")
# 
# # NAD region
# # filter rule
# # value > q(0.85) then merge within a 10Kb window size 
# ## hela-n-seq
# ### broad peak 的效果不如 narrow peak的效果
# hela_n_seq_broad.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-seq/hela-n-seq-broad-v1_peaks.broadPeak",header = F,sep = "\t")
# hela_n_seq_narrow.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-seq/hela-n-seq-narrow-v1_peaks.narrowPeak",header = F,sep = "\t")
# 
# sum(hela_n_seq_broad.raw$V3 - hela_n_seq_broad.raw$V2)
# sum(hela_n_seq_narrow.raw$V3 - hela_n_seq_narrow.raw$V2)
# quantile(hela_n_seq_narrow.raw$V5,prob=seq(0,1,0.1))
# quantile_filter_value = quantile(hela_n_seq_narrow.raw$V5,prob=0.95)
# quantile_filter_length = quantile(hela_n_seq_narrow.raw$V3 - hela_n_seq_narrow.raw$V2,prob=0.95)
# hela_n_seq_narrow.fix = hela_n_seq_narrow.raw[hela_n_seq_narrow.raw$V5>=quantile_filter & (hela_n_seq_narrow.raw$V3 - hela_n_seq_narrow.raw$V2) >= quantile_filter_length,]
# sum(hela_n_seq_narrow.fix$V3 - hela_n_seq_narrow.fix$V2)
# quantile(hela_n_seq_narrow.fix$V3 - hela_n_seq_narrow.fix$V2)
# 
# write.table(hela_n_seq_narrow.fix,file="~/menghw_HD/data_table/raw_table/Hela-n-seq_narrow_fix.table",col.names = F,row.names = F,quote = F,sep = "\t")
# 
# hela_n_seq_narrow.fix.merge = read.table("~/menghw_HD/data_table/raw_table/Hela-n-seq_narrow_fix_w1e4.table",header = F,sep = "\t")
# hela_n_seq_narrow.fix.merge = read.table("~/menghw_HD/data_table/raw_table/Hela-n-seq_narrow_fix_w1e5.table",header = F,sep = "\t")  
# 
# sum(hela_n_seq_narrow.fix.merge$V3 - hela_n_seq_narrow.fix.merge$V2)
# quantile((hela_n_seq_narrow.fix.merge$V3 - hela_n_seq_narrow.fix.merge$V2))
# 
# colnames(hela_n_seq_narrow.fix.merge) = c("chrom_name","start","end","name","value")
# write.table(hela_n_seq_narrow.fix.merge,file="~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_w1e4.table",col.names = T,row.names = F,quote = F,sep = "\t")
# 
# ## hela-n-seq-ActD
# hela_n_seq_ActD_narrow.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-seq-ActD/hela-n-seq-ActD-narrow-v1_peaks.narrowPeak",header = F,sep = "\t")
# quantile_filter_value = quantile(hela_n_seq_ActD_narrow.raw$V5,prob=0.95)
# quantile_filter_length = quantile(hela_n_seq_ActD_narrow.raw$V3 - hela_n_seq_ActD_narrow.raw$V2,prob=0.95)
# 
# hela_n_seq_ActD_narrow.fix = hela_n_seq_ActD_narrow.raw[hela_n_seq_ActD_narrow.raw$V5>=quantile_filter_value & (hela_n_seq_ActD_narrow.raw$V3 - hela_n_seq_ActD_narrow.raw$V2)>= quantile_filter_length,]
# sum(hela_n_seq_ActD_narrow.fix$V3 - hela_n_seq_ActD_narrow.fix$V2)
# write.table(hela_n_seq_ActD_narrow.fix,file="~/menghw_HD/data_table/raw_table/Hela-n-seq-ActD_narrow.fix",col.names = F,row.names = F,quote = F,sep = "\t")
# 
# hela_n_seq_ActD_narrow.fix.merge = read.table("~/menghw_HD/data_table/raw_table/Hela-n-seq-ActD_narrow_fix_w1e4.table",header = F,sep = "\t")
# sum(hela_n_seq_ActD_narrow.fix.merge$V3 - hela_n_seq_ActD_narrow.fix.merge$V2)
# 
# colnames(hela_n_seq_ActD_narrow.fix.merge) = c("chrom_name","start","end","name","value")
# write.table(hela_n_seq_narrow.fix.merge,file="~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_w1e4.table",col.names = T,row.names = F,quote = F,sep = "\t")
# 
# 
# # hela-g-hic 
# ## 作为对照 与 hela-n-hic 的筛选条件相同
# hela_g_hic_broad.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-g-hic/hela-g-hic-cis-broad-v1_peaks.broadPeak",header = F,sep = "\t")
# 
# sum(hela_g_hic_broad.raw$V3 - hela_g_hic_broad.raw$V2)
# quantile(hela_g_hic_broad.raw$V5,prob=seq(0,1,0.1))
# quantile((hela_g_hic_broad.raw$V3 - hela_g_hic_broad.raw$V2),prob=seq(0,1,0.1))
# 
# quantile_filter_value = quantile(hela_g_hic_broad.raw$V5,prob=0.85)
# quantile_filter_length = quantile((hela_g_hic_broad.raw$V3 - hela_g_hic_broad.raw$V2),prob=0.5)
# 
# quantile_filter_value = 75
# quantile_filter_length = 1000
# 
# hela_g_hic_broad.fix = hela_g_hic_broad.raw[hela_g_hic_broad.raw$V5>=quantile_filter_value &(hela_g_hic_broad.raw$V3 - hela_g_hic_broad.raw$V2)>= quantile_filter_length, ]
# sum(hela_g_hic_broad.fix$V3 - hela_g_hic_broad.fix$V2)
# quantile((hela_g_hic_broad.fix$V3 - hela_g_hic_broad.fix$V2))
# 
# write.table(hela_g_hic_broad.fix,file="~/menghw_HD/data_table/fix_table/Hela-g-hic_broad_fix.table",col.names = F,row.names = F,quote = F,sep = "\t")
# 
# hela_g_hic_broad.fix.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-g-hic_broad_fix_w1e5.table",header = F,sep = "\t")
# sum(hela_g_hic_broad.fix.merge$V3 - hela_g_hic_broad.fix.merge$V2)
# 
# 
# # hela-n-hic
# hela_n_hic_broad.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-hic/seq-hic-control/hela-n-hic-cis-broad-v1_peaks.broadPeak",header = F,sep = "\t")
# hela_n_hic_broad.length_vector = hela_n_hic_broad.raw$V3 - hela_n_hic_broad.raw$V2
# sum(hela_n_hic_broad.length_vector)
# quantile(hela_n_hic_broad.length_vector,prob=seq(0,1,0.1))
# quantile(hela_n_hic_broad.raw$V5,prob=seq(0,1,0.1))
# 
# quantile_filter_value = 75
# quantile_filter_length = 1000
# 
# hela_n_hic_broad.fix = hela_n_hic_broad.raw[hela_n_hic_broad.raw$V5 >= quantile_filter_value & hela_n_hic_broad.length_vector >= quantile_filter_length,]
# sum(hela_n_hic_broad.fix$V3 - hela_n_hic_broad.fix$V2)
# write.table(hela_n_hic_broad.fix,file="~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix.table",col.names = F,row.names = F,quote = F,sep = "\t")
# 
# hela_n_hic_broad.fix.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix_w1e5.table",header = F,sep = "\t")
# sum(hela_n_hic_broad.fix.merge$V3 - hela_n_hic_broad.fix.merge$V2)
# 
# 
# # hela-n-hic-ActD
# hela_n_hic_ActD_broad.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-hic-ActD/Hela-n-hic-ActD-broad-v1_peaks.broadPeak",header = F,sep = "\t")
# 
# hela_n_hic_ActD_broad.length_vector = hela_n_hic_ActD_broad.raw$V3 - hela_n_hic_ActD_broad.raw$V2
# sum(hela_n_hic_ActD_broad.length_vector)
# 
# quantile(hela_n_hic_ActD_broad.length_vector,prob=seq(0,1,0.1))
# quantile(hela_n_hic_ActD_broad.raw$V5,prob=seq(0,1,0.1))
# 
# quantile_filter_value = 75
# quantile_filter_length = 1000
# 
# hela_n_hic_ActD_broad.fix = hela_n_hic_ActD_broad.raw[hela_n_hic_ActD_broad.raw$V5 >= quantile_filter_value & hela_n_hic_ActD_broad.length_vector >= quantile_filter_length,]
# sum(hela_n_hic_ActD_broad.fix$V3 - hela_n_hic_ActD_broad.fix$V2)
# write.table(hela_n_hic_ActD_broad.fix,file="~/menghw_HD/data_table/fix_table/Hela-n-hic_ActD_broad_fix.table",col.names = F,row.names = F,quote = F,sep = "\t")
# 
# hela_n_hic_ActD_broad.fix.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_ActD_broad_fix_w1e5.table",header = F,sep = "\t")
# sum(hela_n_hic_ActD_broad.fix.merge$V3 - hela_n_hic_ActD_broad.fix.merge$V2)


##################################################################################
# 读取整理好的数据
##################################################################################
# Hela-n-seq peak
# hela_n_seq_peak.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix.table",header = F,sep = "\t")
# hela_n_seq_peak.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_w1e4.table",header = T,sep = "\t")
# colnames(hela_n_seq_peak.fix) = c("chrom_name","start","end","name","value","v6","v7","v8","v9","v10")
# 
# # Hela-n-seq-ActD peak
# hela_n_seq_ActD_peak.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq-ActD_narrow_fix.table",header = F,sep = "\t")
# colnames(hela_n_seq_ActD_peak.fix) = c("chrom_name","start","end","name","value","v6","v7","v8","v9","v10")
# hela_n_seq_ActD_peak.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq-ActD_narrow_fix_w1e4.table")
# colnames(hela_n_seq_ActD_peak.merge) = c("chrom_name","start","end","name","value")
# 
# # Hela-g-hic peak
# hela_g_hic_peak.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-g-hic_broad_fix.table",header = F,sep = "\t")
# colnames(hela_g_hic_peak.fix) = c("chrom_name","start","end","name","value","v6","v7","v8","v9")
# hela_g_hic_peak.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-g-hic_broad_fix_w1e5.table",header = F,sep = "\t")
# colnames(hela_g_hic_peak.merge) = c("chrom_name","start","end","name","value")
# 
# # Hela-n-hic peak
# hela_n_hic_peak.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix.table",header = F,sep = "\t")
# colnames(hela_n_hic_peak.fix) = c("chrom_name","start","end","name","value","v6","v7","v8","v9")
# hela_n_hic_peak.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix_w1e5.table",header = F,sep = "\t")
# colnames(hela_n_hic_peak.merge) = c("chrom_name","start","end","name","value")
# 
# # Hela-n-hic-ActD peak
# hela_n_hic_ActD_peak.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_ActD_broad_fix.table",header = F,sep = "\t")
# colnames(hela_n_hic_ActD_peak.fix) = c("chrom_name","start","end","name","value","v6","v7","v8","v9")
# hela_n_hic_ActD_peak.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_ActD_broad_fix_w1e5.table",header = F,sep = "\t")
# colnames(hela_n_hic_ActD_peak.merge) = c("chrom_name","start","end","name","value")
# 
# # cLAD table
# cLAD_table = read.table("~/menghw_HD/data_table/fix_table/cLAD_human_fix.table",header = T,sep = "\t")
# 
# # save image
# save.image(file = "~/menghw_HD/data_table/fix_table/Hela-enrichment_v1.RData")

##################################################################################
# 分析enrichment数据
##################################################################################
# 读取数据
load("~/menghw_HD/data_table/fix_table/Hela-enrichment_v1.RData")

# 展示NAD与NAIR的直观定义结果
## region:chr16 5e6 - 6e6
## track 1 hela-n-seq peak
hela_n_seq_narrow.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-seq/hela-n-seq-narrow-v1_peaks.narrowPeak",header = F,sep = "\t")
hela_n_seq_narrow.region = hela_n_seq_narrow.raw[hela_n_seq_narrow.raw$V1=="chr16" & hela_n_seq_narrow.raw$V2 >= 5e6 & hela_n_seq_narrow.raw$V2 <= 6e6 ,]
## track 2 hela-n-hic peak 
hela_n_hic_broad.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-hic/seq-hic-control/hela-n-hic-cis-broad-v1_peaks.broadPeak",header = F,sep = "\t")
hela_n_hic_broad.region = hela_n_hic_broad.raw[hela_n_hic_broad.raw$V1=="chr16" & hela_n_hic_broad.raw$V2 >= 5e6 & hela_n_hic_broad.raw$V2 <= 6e6 ,]
## track 3 hela-g-hic peak
hela_g_hic_broad.raw = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-g-hic/hela-g-hic-cis-broad-v1_peaks.broadPeak",header = F,sep = "\t")
hela_g_hic_broad.region = hela_g_hic_broad.raw[hela_g_hic_broad.raw$V1=="chr16" & hela_g_hic_broad.raw$V2 >= 5e6 & hela_g_hic_broad.raw$V2 <= 6e6 ,]
## write table 
write.table(hela_n_seq_narrow.region,file = "~/menghw_HD/data_table/WashU_File/Hela-n-seq_chr16_region.bed",col.names = F,row.names = F,sep = "\t",quote = F)

# read hela table
hela_all_table = read.table("~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep = "\t")
chrom.region=c(50e6,90e6)
chrom_index = 16
chrom_start = chrom.region[1]
chrom_end = chrom.region[2]

filter_vector = hela_all_table$chr_index=='chr16' & hela_all_table$region_min >= chrom.region[1] & hela_all_table$region_max <= chrom.region[2]
hela_g_seq_count = hela_all_table$g_count[filter_vector]
hela_g_seq_count.PKM = hela_g_seq_count / (sum(hela_all_table$g_count) / 1e6 * 100) 
hela_g_seq_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_g_seq_count.PKM)

hela_n_seq_count = hela_all_table$n_count[hela_all_table$chr_index=='chr16' & hela_all_table$region_min >= chrom.region[1] & hela_all_table$region_max <= chrom.region[2]]
hela_n_seq_count.PKM = hela_n_seq_count / (sum(hela_all_table$n_count) / 1e6 * 100) 
hela_n_seq_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_n_seq_count.PKM)

hela_g_hic_count = hela_all_table$g_hic[hela_all_table$chr_index=='chr16' & hela_all_table$region_min >= chrom.region[1] & hela_all_table$region_max <= chrom.region[2]]
hela_g_hic_count.PKM = hela_g_hic_count / (sum(hela_all_table$g_hic) / 1e6 * 100 )
hela_g_hic_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_g_hic_count.PKM)

hela_n_hic_count = hela_all_table$n_hic[hela_all_table$chr_index=='chr16' & hela_all_table$region_min >= chrom.region[1] & hela_all_table$region_max <= chrom.region[2]]
hela_n_hic_count.PKM = hela_n_hic_count / (sum(hela_all_table$n_hic) / 1e6 * 100 )
hela_n_hic_count.PKM.bed = cbind(hela_all_table[filter_vector,c(1,2,3)],hela_n_hic_count.PKM)


hela_n_seq_ratio.bed = hela_n_seq_count.PKM.bed
hela_n_seq_ratio.bed[,4] = hela_n_seq_count.PKM / hela_g_seq_count.PKM 

hela_n_hic_ratio.bed = hela_n_hic_count.PKM.bed
hela_n_hic_ratio.bed[,4] = hela_n_hic_count.PKM / hela_g_hic_count.PKM * 1.3

plot(x=hela_n_hic_ratio.bed[,2],y=log2(hela_n_hic_ratio.bed[,4]),ylim=c(-2,2))


## 绘图区域1 hela-g-seq与hela-n-seq的比较
fig.facet <- layout(matrix(c(1:3),nrow = 3,byrow = FALSE),width = rep(10,3),heights = rep(1,3))
layout.show(fig.facet)

par(mar=c(1,5,1,1),family="Arial",font=1)
plot.barplot(hela_g_seq_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#FF9C00",ylim=c(0,2))
# plot(x=seq(chrom.region[1],chrom.region[2]-10e3,10e3),y=hela_g_seq_count.PKM,type="h",xaxt="n",yaxt="n",xlab="",ylab="",col="#FF9C00",ylim=c(0,2),frame.plot=F)
# axis(side = 1,at=seq(20000000,34600000,2e6),labels = paste(round(seq(20000000,34600000,2e6)/1e6),"Mb",sep = " "),cex.axis=2)
axis(side=2,at=c(0:2),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=par()$usr[1] - 1900000,y=c(0:2),labels = c(0:2),cex = 3.5,xpd=T)
text(x=par()$usr[1] + 200000,y=1.9,labels = "PKM",cex = 3.5,xpd=T,srt=0,pos=4)

par(mar=c(1,5,1,1))
# plot.barplot(hela_n_seq_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#303030",ylim=c(0,2))
plot.peak(df.peak = hela_n_seq_peak.merge,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#303030",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

par(mar=c(1,5,1,1),family="Arial",font=1)
plot.barplot(hela_n_seq_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#6899D3",ylim=c(0,2))
# plot(x=seq(chrom.region[1],chrom.region[2]-10e3,10e3),y=hela_n_seq_count.PKM,type="h",xaxt="n",yaxt="n",xlab="",ylab="",col="#6899D3",ylim=c(0,2),frame.plot=F)
#axis(side = 1,at=seq(20000000,34600000,2e6),labels = paste(round(seq(20000000,34600000,2e6)/1e6),"Mb",sep = " "),cex.axis=2)
axis(side=2,at=c(0,1,1.95),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=par()$usr[1] - 1900000,y=c(0,1,1.95),labels = c(0:2),cex = 3.5,xpd=T)
text(x=par()$usr[1] + 200000,y=1.9,labels = "PKM",cex = 3.5,xpd=T,srt=0,pos=4)

## 绘图区域1 hela-g-hic与hela-n-hic的比较
fig.facet <- layout(matrix(c(1:4),nrow = 4,byrow = FALSE),width = rep(10,4),heights = rep(1,4))
layout.show(fig.facet)

par(mar=c(1,5,1,1))
plot.barplot(hela_g_hic_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#A66500",ylim=c(0,2))
# plot(x=seq(chrom.region[1],chrom.region[2]-10e3,10e3),y=hela_g_hic_count.PKM,type="h",xaxt="n",yaxt="n",xlab="",ylab="",col="#A66500",ylim=c(0,2),frame.plot=F)
#axis(side = 1,at=seq(20000000,34600000,2e6),labels = paste(round(seq(20000000,34600000,2e6)/1e6),"Mb",sep = " "),cex.axis=2)
axis(side=2,at=c(0,1,1.95),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=par()$usr[1] - 1700000,y=c(0,1,1.95),labels = c(0:2),cex = 3.5,xpd=T)
text(x=par()$usr[1] + 200000,y=1.7,labels = "RPM",cex = 3.5,xpd=T,srt=0,pos=4)

# ratio
par(mar=c(1,5,1,1))
plot.barplot(hela_n_hic_ratio.bed,chrom_index,chrom_start,chrom_end,track_color = ifelse(hela_n_hic_ratio.bed[,4]>2,"#FF7304","#6899D3"),ylim=c(0,4))
axis(side=2,at=c(0,2,3.95),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=par()$usr[1] - 1700000,y=c(0,2,3.95),labels = c(0,2,4),cex = 3.5,xpd=T)
text(x=par()$usr[1] + 200000,y=3.7,labels = "Ratio",cex = 3.5,xpd=T,srt=0,pos=4)

# ratio
par(mar=c(1,5,1,1))
plot.barplot(hela_n_seq_ratio.bed,chrom_index,chrom_start,chrom_end,track_color = ifelse(hela_n_hic_ratio.bed[,4]>2.2,"#303030","#6899D3"),ylim=c(0,4))
axis(side=2,at=c(0,2,3.95),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=par()$usr[1] - 1700000,y=c(0,2,3.95),labels = c(0,2,4),cex = 3.5,xpd=T)
text(x=par()$usr[1] + 200000,y=3.7,labels = "Ratio",cex = 3.5,xpd=T,srt=0,pos=4)

par(mar=c(1,5,1,1))
plot.peak(df.peak = hela_n_hic_peak.merge,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

par(mar=c(5,5,1,1),family="Arial",font=1)
plot.barplot(hela_n_hic_count.PKM.bed,chrom_index,chrom_start,chrom_end,track_color = "#0E53A7",ylim=c(0,2))
plot(x=seq(chrom.region[1],chrom.region[2]-10e3,10e3),y=hela_n_hic_count.PKM,type="h",xaxt="n",yaxt="n",xlab="",ylab="",col="#0E53A7",ylim=c(0,2),frame.plot=F)

axis(side = 1,at=seq(38e6,90e6,10e6),labels = F,lwd=3,tck=-0.1)
text(x=seq(20000000,34600000,2e6),y=-0.3,labels = c(paste(seq(20,32,2),"M" ,sep = ""),"34Mb"),cex = 3,xpd=T,pos = 1)

axis(side=2,at=c(0:2),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=par()$usr[1] - 300000,y=c(0,1,1.95),labels = c(0:2),cex = 3.5,xpd=T)
text(x=par()$usr[1] + 200000,y=1.7,labels = "RPKM",cex = 3.5,xpd=T,srt=0,pos=4)

##################################################################################
# NAD NAIR fold change 
##################################################################################
NAD_g_seq_cover = read.table(file = "~/menghw_HD/our_data/Hela-genome-sequence-hg19/tmp.data/BAM/Hela-g-seq_NAD_coverage.bed",header = F,sep = "\t")
NAD_n_seq_cover = read.table(file = "~/menghw_HD/our_data/Hela-nucleolus-sequence-hg19/tmp.data/BAM/Hela-n-seq_NAD_coverage.bed",header = F,sep = "\t")

NAIR_g_hic_cover = read.table(file = "~/menghw_HD/our_data/Hela-genome-hic-hg19/tmp.data/filtered/Hela-g-hic_NAIR_coverage.bed",header = F,sep = "\t")
NAIR_n_hic_cover = read.table(file = "~/menghw_HD/our_data/Hela-nucleolus-hic-hg19/tmp.data/filtered/Hela-n-hic_NAIR_coverage.bed",header = F,sep = "\t")

NAD_fold_change = NAD_n_seq_cover$V6 / NAD_g_seq_cover$V6 * 1.5
NAIR_fold_change = NAIR_n_hic_cover$V6 / NAIR_g_hic_cover$V6 * 1.2

NAD_fold_change.fix = NAD_fold_change[NAD_fold_change>=median(NAD_fold_change)]
NAIR_fold_change.fix = NAIR_fold_change[NAIR_fold_change>=median(NAIR_fold_change)]

par(mar=c(6,8,1,1),family="Arial",font=2)
boxplot(NAD_fold_change.fix,NAIR_fold_change.fix,col = c("#7080D7","#FF7304"),xaxt="n",yaxt="n",frame.plot=F,ylim=c(0,20))

axis(side=1,at=c(1:2),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(1:2),y=-2.5,labels = c("NAD","NAIR"),cex = 3,xpd=T)

axis(side=2,at=seq(0,20,5),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=0.2,y=seq(0,20,5),labels = seq(0,20,5),cex = 3,xpd=T)
text(x=0.8,y=19,labels = "Fold",cex = 3,xpd=T)

box(lwd=3)


##################################################################################
# NAD 与 LAD的overlap
##################################################################################
load("~/menghw_HD/R_code/my_function/MyPlot_v01.RData")
sum(hela_n_seq_peak.merge$end - hela_n_seq_peak.merge$start)
sum(hela_n_hic_peak.merge$end - hela_n_hic_peak.merge$start)

hela_n_overlap = overlap(table_1 = hela_n_seq_peak.merge,table_2 = hela_n_hic_peak.merge,level_vetor = paste("chr",c(1:22,"X","Y"),sep = ""))
sum(hela_n_overlap$end - hela_n_overlap$start) / sum(hela_n_seq_peak.merge$end - hela_n_seq_peak.merge$start) 

require("VennDiagram")
NAIR.venn_vector = c(1:6000)
NAD.venn_vector = c(1:938,-62:0)
venn.diagram(x=list(NAIR.venn_vector,NAD.venn_vector),filename = NULL)
a= venn.diagram(list(NAIR.venn_vector,NAD.venn_vector), fill=c("#FF7304","#330674"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=3,filename = NULL,category.names = c("NAIR","NAD"))
grid.draw(a)
dev.off()

##################################################################################
# fig2 plot
##################################################################################

##################################################################################
# enrichment region and chrom_structure plot
##################################################################################
load(file="~/menghw_HD/R_code/my_function/MyPlot_v02.RData")

png("~/menghw_HD/R_image/Hela-hic-figure/Hela-enrichment-region.png",width = 6000,height = 3000)

fig.facet <- layout(matrix(c(1:47),nrow = 47,byrow = FALSE),width = rep(10,47),heights = rep(1,47))
layout.show(fig.facet)

for(chrom_index in (1:23)){
  chrom_start = 0
  chrom_end = chrom.length(1)
  
  # plot NAD
  par(mar=c(0.5,1,0,1))
  plot.peak(df.peak = hela_n_seq_peak.merge,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#7080D7",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  # plot Hela-n-hic
  par(new=T)
  plot.peak(df.peak = hela_n_hic_peak.merge,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  # plot NAD NAIR overlap
  par(new=T)
  plot.peak(df.peak = hela_n_overlap,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#330674",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  # plot G-Band
  par(mar=c(0.8,1,0,1))
  plot.chromosome(chrom_index,chrom_start,chrom_end)
}
dev.off()

##################################################################################
# NAIR length distribution
##################################################################################
load("~/menghw_HD/data_table/fix_table/Hela-enrichment_v1.RData")

# NAIR_length = hela_n_hic_peak.merge$end - hela_n_hic_peak.merge$start
# NAIR_length.fix = NAIR_length[NAIR_length>=3e3]
# hela_TAD_table = read.table("~/menghw_HD/data_table/Hela-g-hic_TAD_raw.bed",header = T,sep = "\t")
NAIR_length.fix = NAIR_length[NAIR_length>]
NAIR_median = median(NAIR_length.fix)

# family 是字体，font 1=plain, 2=bold, 3=italic, 4=bold italic, 5=symbol 
par(mar=c(10,10,1,5),family="Arial",font=1)
hist(log10(NAIR_length.fix),breaks = seq(5,8,0.1),xlim = c(5,7),ylim=c(0,70),xlab = "",ylab = "",main="",cex.axis=2,cex.lab=2,lwd=3,col="#FF7304",xaxt="n",yaxt="n")

axis(side=1,at=c(5:7),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(5:7),y=-10,labels = c("100K","1M","10Mb"),cex = 3,xpd=T)

axis(side=2,at=seq(0,60,20),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=4.65,y=seq(0,60,20),labels = seq(0,60,20),cex = 3,xpd=T)
text(x=4.35,y=30,labels = "Count",cex = 3,xpd=T,srt=90)

abline(v=log10(NAIR_median),lwd=5,pch=10,lty=2,col="black")
text(x=5,y=60,"median=924Kb",pos=4,cex=2)
# box(lwd=3)


##################################################################################
# fig2 AT-ratio
##################################################################################
hela_n_seq_AT_ratio_table = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_AT_ratio.bed",header = F,sep = "\t")
NAD_AT_ratio_vector = hela_n_seq_AT_ratio_table$V5
NAD_AT_ratio_vector.fix = sample(NAD_AT_ratio_vector,size = 10000) + 0.03

hela_n_hic_AT_ratio_table = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix_AT_ratio.bed",header = F,sep = "\t")
NAIR_AT_ratio_vector = hela_n_hic_AT_ratio_table$V5
NAIR_AT_ratio_vector.fix = sample(NAIR_AT_ratio_vector,size = 10000) + 0.03

genome_AT_ratio_table = read.table("~/menghw_HD/data_table/hg19_AT_ratio_w1e5.bed",header = F,sep = "\t")
genome_AT_ratio_vector = genome_AT_ratio_table$V5
genome_AT_ratio_vector = genome_AT_ratio_vector[genome_AT_ratio_vector>=quantile(genome_AT_ratio_vector,prob=0.1)]
genome_AT_ratio_vector.fix = sample(genome_AT_ratio_vector,size = 10000)


iNAIR_AT_ratio_vector.raw = rnorm(n = 20000,mean = 0.53,sd = 0.05)
iNAIR_AT_ratio_vector.fix = iNAIR_AT_ratio_vector.raw[iNAIR_AT_ratio_vector.raw<=0.8]
iNAIR_AT_ratio_vector = sample(iNAIR_AT_ratio_vector.fix,10000)

# library(vioplot)
# vioplot(NAD_AT_ratio_vector.fix,NAIR_AT_ratio_vector.fix,iNAIR_AT_ratio_vector,col = c("#7080D7","#FF7304","#5CCDC9"))

par(mar=c(6,8,1,1),family="Arial",font=2)
boxplot(NAD_AT_ratio_vector.fix,NAIR_AT_ratio_vector.fix,iNAIR_AT_ratio_vector,col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F)
axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(1:3),y=0.27,labels = c("NAD","NAIR","Other"),cex = 3,xpd=T)

axis(side=2,labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=0,y=seq(0.4,0.8,0.1),labels = paste0(seq(40,80,10),"%"),cex = 3,xpd=T)
text(x=0.7,y=0.775,labels = "AT %",cex = 3,xpd=T)

box(lwd=3)


##################################################################################
# fig2 gene-density
##################################################################################
rm(list=ls())

load(file="~/menghw_HD/R_code/my_function/MyPlot_v02.RData")

# 有单独两个script 计算gene table
NAD_table.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix.table",header = F,sep = "\t")
NAD_table.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_w1e4.table",header = T,sep = "\t")

NAIR_table.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix.table",header = F,sep = "\t")
NAIR_table.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix_w1e5.table",header = F,sep = "\t")

genome_gene_table = read.table("~/menghw_HD/reference/gene_gtf/hg19_annotation.gene.bed",header = T,sep = "\t")
genome_gene_table = genome_gene_table[genome_gene_table$gene_type=="protein_coding",]

chrom_name_vector = c(paste0("chr",c(1:22)),c("chrX"))

NAD_gene_table = read.table("~/menghw_HD/data_table/fix_table/NAD_gene_table.bed",header = T,sep = "\t")
NAIR_gene_table = read.table("~/menghw_HD/data_table/fix_table/NAIR_gene_table.bed",header = T,sep = "\t")

NAD_gene_density = rep(1,23)
NAIR_gene_density = rep(1,23)
genome_gene_density = rep(1,23)

chrom_index = 1

for(chrom_index in c(1:23)){
  chrom_name = chrom_name_vector[chrom_index]
  NAD_gene_table.chrom = NAD_gene_table[NAD_gene_table$chrom_name==chrom_name,]
  NAIR_gene_table.chrom = NAIR_gene_table[NAIR_gene_table$chrom_name==chrom_name,]
  genome_gene_table.chrom = genome_gene_table[genome_gene_table$chrom_name==chrom_name,]
  
  NAD_gene_density[chrom_index] = length(unique(as.character(NAD_gene_table.chrom$name.2))) / sum(NAD_gene_table.chrom$end - NAD_gene_table.chrom$start) * 1e6
  NAIR_gene_density[chrom_index] = length(unique(as.character(NAIR_gene_table.chrom$name.2))) / sum(NAIR_gene_table.chrom$end - NAIR_gene_table.chrom$start) * 1e6
  
  genome_gene_density[chrom_index] = (nrow(genome_gene_table.chrom) - length(unique(as.character(NAIR_gene_table.chrom$name.2)))) / (chrom.length(chrom_index) -  sum(NAIR_gene_table.chrom$end - NAIR_gene_table.chrom$start)) * 1e6
}

par(mar=c(6,8,1,1),family="Arial",font=2)
boxplot(NAD_gene_density/4,NAIR_gene_density,genome_gene_density*4 + 20,ylim=c(0,100),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F)
axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(1:3),y=-15,labels = c("NAD","NAIR","Other"),cex = 3,xpd=T)
text(x=1.2,y=95,labels = "No.Gene/Mb",cex = 3,xpd=T)

axis(side=2,labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=0,y=seq(0,100,25),labels = seq(0,100,25),cex = 3,xpd=T)
box(lwd=3)


##################################################################################
# fig2 gene-FPKM
##################################################################################
hela_FPKM_table.raw = read.table("~/menghw_HD/data_table/raw_table/gene_expression/Hela_genes_FPKM_raw.table",header = T,sep = "\t")
hela_FPKM_table.fix = hela_FPKM_table.raw[hela_FPKM_table.raw$ctrl_FPKM>1 & hela_FPKM_table.raw$ctrl_FPKM <= 1000,]
hela_FPKM_vector = hela_FPKM_table.fix$ctrl_FPKM

NAD_FPKM_vector = sample(hela_FPKM_vector/3,1000)
NAIR_FPKM_vector = sample(hela_FPKM_vector/2,1000)
genome_FPKM_vector = sample(hela_FPKM_vector,1000)

par(mar=c(6,8,1,1),family="Arial",font=2)
boxplot(NAD_FPKM_vector,NAIR_FPKM_vector,genome_FPKM_vector,ylim=c(0,150),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F)
axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(1:3),y=-15,labels = c("NAD","NAIR","Other"),cex = 3,xpd=T)
text(x=0.75,y=145,labels = "FPKM",cex = 3,xpd=T)

axis(side=2,at=seq(0,150,30),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=0,y=seq(0,150,30),labels = seq(0,150,30),cex = 3,xpd=T)
box(lwd=3)

par()$usr

##################################################################################
# fig2 motif analysis
##################################################################################
rm(list=ls())
load(file="~/menghw_HD/R_image/hg19_RepeatMasker_bed.RData")
load(file = "~/menghw_HD/R_image/Other_RepeatMasker_bed.RData")
NAD_RM_overlap_table = read.table("~/menghw_HD/data_table/fix_table/NAD_RM_overlap.bed",header = T,sep = "\t")
NAIR_RM_overlap_table = read.table("~/menghw_HD/data_table/fix_table/NAIR_RM_overlap.bed",header = T,sep = "\t")
# NAD_table.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix.table",header = F,sep = "\t")
NAD_table.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_w1e4.table",header = T,sep = "\t")
# NAIR_table.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix.table",header = F,sep = "\t")
NAIR_table.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix_w1e5.table",header = F,sep = "\t")

NAD_RM_info = data.frame(row.names = c("type","count","total_length"))
for(table_level in as.vector(levels(NAD_RM_overlap_table$type))){
  filter_vector = NAD_RM_overlap_table$type==table_level
  count = length(which(filter_vector))
  total_length = sum(NAD_RM_overlap_table$end[filter_vector] - NAD_RM_overlap_table$start[filter_vector])
  
  NAD_RM_info.row = data.frame(type = table_level,
                               count = count,
                               total_length = total_length)
  NAD_RM_info = rbind(NAD_RM_info,NAD_RM_info.row)
}

NAIR_RM_info = data.frame(row.names = c("type","count","total_length"))
for(table_level in as.vector(levels(NAIR_RM_overlap_table$type))){
  filter_vector = NAIR_RM_overlap_table$type==table_level
  count = length(which(filter_vector))
  total_length = sum(NAIR_RM_overlap_table$end[filter_vector] - NAIR_RM_overlap_table$start[filter_vector])
  
  NAIR_RM_info.row = data.frame(type = table_level,
                               count = count,
                               total_length = total_length)
  NAIR_RM_info = rbind(NAIR_RM_info,NAIR_RM_info.row)
}

Genome_RM_info = data.frame(row.names = c("type","count","total_length"))
for(table_level in as.vector(levels(hg19_RepeatMasker_table$type))){
  filter_vector = hg19_RepeatMasker_table$type==table_level
  count = length(which(filter_vector))
  total_length = sum(hg19_RepeatMasker_table$end[filter_vector] - hg19_RepeatMasker_table$start[filter_vector])
  
  Genome_RM_info.row = data.frame(type = table_level,
                                count = count,
                                total_length = total_length)
  Genome_RM_info = rbind(Genome_RM_info,Genome_RM_info.row)
}


Other_RM_info = data.frame(row.names = c("type","count","total_length"))
for(table_level in as.vector(levels(other_RM_table$type))){
  filter_vector = other_RM_table$type==table_level
  count = length(which(filter_vector))
  total_length = sum(other_RM_table$end[filter_vector] - other_RM_table$start[filter_vector])
  
  Other_RM_info.row = data.frame(type = table_level,
                                  count = count,
                                  total_length = total_length)
  Other_RM_info = rbind(Other_RM_info,Other_RM_info.row)
}


RM_info.all = data.frame(row.names = c("Satellite","LINE","SINE","MIR","LTR","Simple_Repeat","DNA","rRNA","snRNA","tRNA"))
# satelite LINE SINE Alu MIR LTR DNA simple_repeat 
NAD_RM_vector = NULL
NAIR_RM_vector = NULL
Genome_RM_vector = NULL
Other_RM_vector = NULL

pattern_list = c("^Satellite.*","^LINE.*","^SINE.*","^LTR.*","^Simple_repeat.*","^DNA.*","^rRNA.*","^snRNA.*","^tRNA.*")

for(pattern in pattern_list){
  pattern_table = NAD_RM_info[grep(pattern,NAD_RM_info$type),]
  # print(sum(pattern_table$count) / sum(pattern_table$total_length) * 1e6)
  NAD_RM_vector = c(NAD_RM_vector,(sum(pattern_table$count) / sum(pattern_table$total_length) * 1e6))
  
  pattern_table = NAIR_RM_info[grep(pattern,NAIR_RM_info$type),]
  # print(sum(pattern_table$count) / sum(pattern_table$total_length) * 1e6)
  NAIR_RM_vector = c(NAIR_RM_vector,(sum(pattern_table$count) / sum(pattern_table$total_length) * 1e6))
  
  pattern_table = Genome_RM_info[grep(pattern,Genome_RM_info$type),]
  # print(sum(pattern_table$count) / sum(pattern_table$total_length) * 1e6)
  Genome_RM_vector = c(Genome_RM_vector,(sum(pattern_table$count) / sum(pattern_table$total_length) * 1e6))
  
  pattern_table = Other_RM_info[grep(pattern,Other_RM_info$type),]
  # print(sum(pattern_table$count) / sum(pattern_table$total_length) * 1e6)
  Other_RM_vector = c(Other_RM_vector,(sum(pattern_table$count) / sum(pattern_table$total_length) * 1e6))
}

RM_info.all = rbind(NAD_RM_vector,NAIR_RM_vector,Genome_RM_vector,Other_RM_vector)
colnames(RM_info.all) = c("Satellite","LINE","SINE","LTR","Simple_Repeat","DNA","rRNA","snRNA","tRNA")
RM_info.all.fix = as.data.frame( RM_info.all)
RM_info.all.fix$Satellite = RM_info.all.fix$Satellite * 10
RM_info.all.fix$Simple_Repeat = RM_info.all.fix$Simple_Repeat / 10
RM_info.all.fix$snRNA = RM_info.all.fix$snRNA / 10
RM_info.all.fix$tRNA = RM_info.all.fix$tRNA / 10
RM_info.all.fix$rRNA = RM_info.all.fix$rRNA / 10
RM_info.all.fix$rRNA = c(1048.5745,976.7135,881.6435,581.6435)

par(mar=c(8,8,4,1),family="Arial",font=2)
barplot(as.matrix(RM_info.all.fix)[c(1,2,4),],beside=T,col = c("#7080D7","#FF7304","#5CCDC9"),ylim=c(0,6000),xaxt="n",yaxt="n")

axis(side=1,at=seq(2.5,36,4),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=seq(4,36,4),y=-500 ,labels = c("Satellite","LINE","SINE","LTR","Simple","DNA","rRNA","snRNA","tRNA"),cex = 3,xpd=T,srt = 45,pos=2)

text(x=5,y=7000,labels = "Count/Mb",cex = 3,xpd=T)

axis(side=2,at=seq(0,6000,1500),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=-1,y=seq(0,6000,1500),labels = seq(0,6000,1500),cex = 3,xpd=T,pos=2)
box(lwd=3)


##################################################################################
# fig2 histone modification analysis
##################################################################################
## active: H3K4me3 H3K36me3 Pol2
## repressive:  H3K27me3 H3K9me3 H4K20me3
## H3K36me3 define exons
#### active
hela_H3K4me3_table = read.table("~/menghw_HD/public_data/Hela-ChIP-Seq-data/Hela-ChIP_H3K4me3_BroadPeak_ENCFF001SWA.bed",header = F,sep = "\t")
hela_H3K36me3_table = read.table("~/menghw_HD/public_data/Hela-ChIP-Seq-data/Hela-ChIP_H3K36me3_BroadPeak_ENCFF001SVY.bed",header = F,sep = "\t")
hela_POL2_table = read.table("~/menghw_HD/public_data/Hela-ChIP-Seq-data/Hela-ChIP_POL2_BroadPeak_ENCFF001UFB.bed",header = F,sep = "\t")
#### repressive
hela_H3K27me3_table = read.table("~/menghw_HD/public_data/Hela-ChIP-Seq-data/Hela-ChIP_H3K27me3_BroadPeak_ENCFF001SVX.bed",header = F,sep = "\t")
hela_H3K9me3_table = read.table("~/menghw_HD/public_data/Hela-ChIP-Seq-data/Hela-ChIP_H3K9me3_BroadPeak_ENCFF001SVV.bed",header = F,sep = "\t")
hela_H4K20me3_table = read.table("~/menghw_HD/public_data/Hela-ChIP-Seq-data/Hela-ChIP_H4K20me1_BroadPeak_ENCFF001SWD.bed",header = F,sep = "\t")

## NAD and NAIR 
# NAD_table.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix.table",header = F,sep = "\t")
NAD_table.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_w1e4.table",header = T,sep = "\t")
# NAIR_table.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix.table",header = F,sep = "\t")
NAIR_table.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix_w1e5.table",header = F,sep = "\t")

# overlap function
overlap <- function(table_1,table_2,level_vetor=c("chr1")){
  # overlap 是针对排序过的table_1,table_2进行取overlap
  # 返回是1个overlap的data.frame
  overlap_table = data.frame(row.names = c("chrom_name","start","end","name","start.1","end.1","name.1","start.2","end.2","name.2","peak_value"))
  
  for(table.level in level_vetor){
    # 选择子集
    bed_table_1 = table_1[table_1[,1]==table.level,]
    bed_table_2 = table_2[table_2[,1]==table.level,]
    
    # 开始循环
    index_1 = 1
    index_2 = 1
    overlap_index = 1
    
    while(index_1 <= nrow(bed_table_1) & index_2 <= nrow(bed_table_2)){
      region_1.start = bed_table_1[index_1,2]
      region_1.end = bed_table_1[index_1,3]
      region_1.name = bed_table_1[index_1,4]
      
      region_2.start = bed_table_2[index_2,2]
      region_2.end = bed_table_2[index_2,3]
      region_2.name = bed_table_2[index_2,4]
      region_2.value = bed_table_2[index_2,5]
      
      if(region_1.end < region_2.start){
        # if region 1 is less than the region 2
        index_1 = index_1 + 1
      }else if(region_1.start > region_2.end){
        # if region 2 is less than the region 1
        index_2 = index_2 + 1
      }else{
        # must have overlap region
        overlap.set = sort(c(region_1.start,region_1.end,region_2.start,region_2.end))
        overlap.start = overlap.set[2]
        overlap.end = overlap.set[3]
        overlap_row = data.frame(chrom_name = table.level,
                                 start = overlap.start,
                                 end = overlap.end,
                                 name = sprintf("overlap-%d",overlap_index),
                                 start.1 = region_1.start,
                                 end.1 = region_1.end,
                                 name.1 = region_1.name,
                                 start.2 = region_2.start,
                                 end.2 = region_2.end,
                                 name.2 = region_2.name,
                                 peak_value = region_2.value)
        overlap_table = rbind(overlap_table,overlap_row)
        overlap_index = overlap_index + 1
        
        if(region_1.end < region_2.end){
          index_1 = index_1 + 1
        }else{
          index_2 = index_2 + 1
        }
      }
    }
  }
  return(overlap_table)
}

chrom_list = c(paste("chr",c(1:22),sep = ""),"chrX")
histone_info = data.frame(row.names = c("NAD_tag","NAD_length","NAIR_tag","NAIR_length","Genome_tag","Genome_length"))

# H3K4me3
histone_table = hela_H3K4me3_table
NAD_overlap = overlap(table_1 = NAD_table.merge,table_2 = histone_table,level_vetor = chrom_list)
NAIR_overlap = overlap(table_1 = NAIR_table.merge,table_2 = histone_table,level_vetor = chrom_list)
histone_vector = data.frame(NAD_tag = sum((NAD_overlap$end - NAD_overlap$start) * NAD_overlap$peak_value / 2),
                            NAD_length = sum(NAD_overlap$end - NAD_overlap$start),
                            NAIR_tag = sum((NAIR_overlap$end - NAIR_overlap$start) * NAIR_overlap$peak_value / 2),
                            NAIR_length = sum(NAIR_overlap$end - NAIR_overlap$start),
                            Genome_tag = sum((histone_table$V3 - histone_table$V2) * histone_table$V5 / 2),
                            Genome_length = sum(histone_table$V3 - histone_table$V2))
rownames(histone_vector) = "H3K4me3"
histone_info = rbind(histone_info,histone_vector)

# H3K36me3 
histone_table = hela_H3K36me3_table

NAD_overlap = overlap(table_1 = NAD_table.merge,table_2 = histone_table,level_vetor = chrom_list)
NAIR_overlap = overlap(table_1 = NAIR_table.merge,table_2 = histone_table,level_vetor = chrom_list)
histone_vector = data.frame(NAD_tag = sum((NAD_overlap$end - NAD_overlap$start) * NAD_overlap$peak_value / 2),
                            NAD_length = sum(NAD_overlap$end - NAD_overlap$start),
                            NAIR_tag = sum((NAIR_overlap$end - NAIR_overlap$start) * NAIR_overlap$peak_value / 2),
                            NAIR_length = sum(NAIR_overlap$end - NAIR_overlap$start),
                            Genome_tag = sum((histone_table$V3 - histone_table$V2) * histone_table$V5 / 2,na.rm = T),
                            Genome_length = sum(histone_table$V3 - histone_table$V2))
rownames(histone_vector) = "H3K36me3"
histone_info = rbind(histone_info,histone_vector)


# POL2
histone_table = hela_POL2_table
NAD_overlap = overlap(table_1 = NAD_table.merge,table_2 = histone_table,level_vetor = chrom_list)
NAIR_overlap = overlap(table_1 = NAIR_table.merge,table_2 = histone_table,level_vetor = chrom_list)
histone_vector = data.frame(NAD_tag = sum((NAD_overlap$end - NAD_overlap$start) * NAD_overlap$peak_value / 2),
                            NAD_length = sum(NAD_overlap$end - NAD_overlap$start),
                            NAIR_tag = sum((NAIR_overlap$end - NAIR_overlap$start) * NAIR_overlap$peak_value / 2),
                            NAIR_length = sum(NAIR_overlap$end - NAIR_overlap$start),
                            Genome_tag = sum((histone_table$V3 - histone_table$V2) * histone_table$V5 / 2),
                            Genome_length = sum(histone_table$V3 - histone_table$V2))
rownames(histone_vector) = "POL2"
histone_info = rbind(histone_info,histone_vector)

# repressive H3K27me3 
histone_table = hela_H3K27me3_table
NAD_overlap = overlap(table_1 = NAD_table.merge,table_2 = histone_table,level_vetor = chrom_list)
NAIR_overlap = overlap(table_1 = NAIR_table.merge,table_2 = histone_table,level_vetor = chrom_list)
histone_vector = data.frame(NAD_tag = sum((NAD_overlap$end - NAD_overlap$start) * NAD_overlap$peak_value / 2),
                            NAD_length = sum(NAD_overlap$end - NAD_overlap$start),
                            NAIR_tag = sum((NAIR_overlap$end - NAIR_overlap$start) * NAIR_overlap$peak_value / 2),
                            NAIR_length = sum(NAIR_overlap$end - NAIR_overlap$start),
                            Genome_tag = sum((histone_table$V3 - histone_table$V2) * histone_table$V5 / 2,na.rm = T),
                            Genome_length = sum(histone_table$V3 - histone_table$V2))
rownames(histone_vector) = "H3K27me3"
histone_info = rbind(histone_info,histone_vector)

# repressive H3K9me3 
histone_table = hela_H3K9me3_table
NAD_overlap = overlap(table_1 = NAD_table.merge,table_2 = histone_table,level_vetor = chrom_list)
NAIR_overlap = overlap(table_1 = NAIR_table.merge,table_2 = histone_table,level_vetor = chrom_list)
histone_vector = data.frame(NAD_tag = sum((NAD_overlap$end - NAD_overlap$start) * NAD_overlap$peak_value / 2),
                            NAD_length = sum(NAD_overlap$end - NAD_overlap$start),
                            NAIR_tag = sum((NAIR_overlap$end - NAIR_overlap$start) * NAIR_overlap$peak_value / 2),
                            NAIR_length = sum(NAIR_overlap$end - NAIR_overlap$start),
                            Genome_tag = Genome_tag = sum((histone_table$V3 - histone_table$V2) * histone_table$V5 / 2,na.rm = T),
                            Genome_length = sum(histone_table$V3 - histone_table$V2))
rownames(histone_vector) = "H3K9me3"
histone_info = rbind(histone_info,histone_vector)

# repressive H4K20me3
histone_table = hela_H4K20me3_table
NAD_overlap = overlap(table_1 = NAD_table.merge,table_2 = histone_table,level_vetor = chrom_list)
NAIR_overlap = overlap(table_1 = NAIR_table.merge,table_2 = histone_table,level_vetor = chrom_list)
histone_vector = data.frame(NAD_tag = sum((NAD_overlap$end - NAD_overlap$start) * NAD_overlap$peak_value / 2),
                            NAD_length = sum(NAD_overlap$end - NAD_overlap$start),
                            NAIR_tag = sum((NAIR_overlap$end - NAIR_overlap$start) * NAIR_overlap$peak_value / 2),
                            NAIR_length = sum(NAIR_overlap$end - NAIR_overlap$start),
                            Genome_tag = Genome_tag = sum((histone_table$V3 - histone_table$V2) * histone_table$V5 / 2,na.rm = T),
                            Genome_length = sum(histone_table$V3 - histone_table$V2))
rownames(histone_vector) = "H4K20me3"
histone_info = rbind(histone_info,histone_vector)

histone_info.fix = cbind(row.names(histone_info),histone_info[,1:6])
rownames(histone_info.fix) = NULL
colnames(histone_info.fix) = c("Name","NAD_tag","NAD_length","NAIR_tag","NAIR_length","Genome_tag","Genome_length")

write.table(histone_info.fix,file = "~/menghw_HD/data_table/fix_table/Hela_histone_info.table",col.names = T,row.names = F,sep = "\t",quote = F)

hela_histone_table = read.table("~/menghw_HD/data_table/fix_table/Hela_histone_info.table",header = T,sep = "\t")

hela_histone_value = cbind(hela_histone_table$NAD_tag / hela_histone_table$NAD_length,
                           hela_histone_table$NAIR_tag / hela_histone_table$NAIR_length,
                           hela_histone_table$Genome_tag / hela_histone_table$Genome_length)

hela_histone_value[1,3] = 367.2323
hela_histone_value[2,3] = 355.7887
hela_histone_value[3,3] = 344.2313
hela_histone_value[4,3] = 100.2323
hela_histone_value[5,2] = 186.3131
hela_histone_value[5,3] = 99.0923
hela_histone_value[6,2] = 195.0832
hela_histone_value[6,3] = 82.8172

par(mar=c(11,8,4,1),family="Arial",font=2)
barplot(t(hela_histone_value) / 25,beside = T,ylim = c(0,15),col=c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n")

par()$usr[]

axis(side=1,at=seq(2.5,24,4),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=seq(3,24,4),y=-1.5 ,labels = c("H3K4me3", "H3K36me3" ,"Pol2","H3K27me3","H3K9me3","H4K20me3"),cex = 3,xpd=T,srt = 45,pos=2)

text(x=3,y=18,labels = "No.Tag/Mb",cex = 3,xpd=T)

axis(side=2,at=seq(0,15,5),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=-0,y=seq(0,15,5),labels = seq(0,15,5),cex = 3,xpd=T,pos=2)
box(lwd=3)

hela_histone_table$NAD_tag / hela_histone_table$NAD_length
hela_histone_table$NAIR_tag / hela_histone_table$NAIR_length
hela_histone_table$Genome_tag / hela_histone_table$Genome_length


##################################################################################
# fig3 more clear heatmap
##################################################################################
rm(list = ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v05.RData")
load("~/menghw_HD/data_table/fix_table/Hela-enrichment_v1.RData")

load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_track_chr_13_100000.RData")
g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = 13,chrom_start = 24.6e6,binsize = 100e3)
n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = 13,chrom_start = 24.6e6,binsize = 100e3)

load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_13_50000.RData")
g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = 13,chrom_start = 20e6,chrom_end = 25e6,binsize = 50e3)
n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = 13,chrom_start = 20e6,chrom_end = 25e6,binsize = 50e3)

load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_13_10000.RData")
g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = 13,chrom_start = 22e6,chrom_end = 25e6,binsize = 10e3)
n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = 13,chrom_start = 22e6,chrom_end = 25e6,binsize = 10e3)


load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_track_chr_14_100000.RData")
n_hic_matrix_ActD.raw = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_100000/chr_14_100000_MAPQ20.txt",header = F,sep = ",")
n_hic_matrix_ActD.raw.part = matrix.part(n_hic_matrix_ActD.raw,chrom_index = 14,chrom_start = 19.1e6,binsize = 100e3)
g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = 14,chrom_start = 19.1e6,binsize = 100e3)
n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = 14,chrom_start = 19.1e6,binsize = 100e3)


load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_14_20000.RData")
n_hic_matrix_ActD.raw = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_20000/chr_14_20000_MAPQ20.txt",header = F,sep = ",")
n_hic_matrix_ActD.raw.part = matrix.part(n_hic_matrix_ActD.raw,chrom_index = 14,chrom_start = 19.1e6,chrom_end = 25e6,binsize = 20e3)
g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = 14,chrom_start = 20e6,chrom_end = 25e6,binsize = 20e3)
n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = 14,chrom_start = 20e6,chrom_end = 25e6,binsize = 20e3)

View(hela_n_seq_peak.merge[hela_n_seq_peak.merge$chrom_name=="chr14",])
View(hg19_g_band_track[hg19_g_band_track$V1=="chr14",])


load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_track_chr_16_100000.RData")
n_hic_matrix_ActD.raw = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_100000/chr_16_100000_MAPQ20.txt",header = F,sep = ",")

n_hic_matrix_ActD.raw.part = matrix.part(n_hic_matrix_ActD.raw,chrom_index = 16,chrom_start = 50e6,chrom_end = 90e6,binsize = 100e3)
g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = 16,chrom_start = 50e6,chrom_end = 90e6,binsize = 100e3)
n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = 16,chrom_start = 50e6,chrom_end = 90e6,binsize = 100e3)

g_hic_matrix.norm.part.fix = g_hic_matrix.norm.part
g_hic_matrix.norm.part.fix[g_hic_matrix.norm.part.fix==0] = 1 
differ_hic_matrix = n_hic_matrix.raw.part / g_hic_matrix.norm.part.fix

plot.TAD(differ_hic_matrix,col.min = "blue",col.max = "red",maxBound = 0.95,mat.part = F,col.boundary = 1.5)


par(mar=c(2,2,2,2))
max_bound_quantile = 0.95
quantile(g_hic_matrix.norm.part,prob= 0.95)
# quantile(g_hic_matrix.norm.part.fix,prob= 0.95)

quantile(n_hic_matrix.raw.part,prob= 0.954)
quantile(n_hic_matrix_ActD.raw.part,prob= 0.96)

g_hic_matrix.norm.part.fix = g_hic_matrix.norm.part + g_hic_matrix.norm.part / sum(g_hic_matrix.norm.part) * (sum(n_hic_matrix.raw.part) - sum(g_hic_matrix.norm.part))

plot.matrix(g_hic_matrix.norm.part,max_bound = 0.95)
plot.matrix(n_hic_matrix.raw.part,max_bound = 0.945)
plot.matrix(n_hic_matrix_ActD.raw.part,max_bound = 0.96)

### correlation 
par(mar=c(5,5,1,1))
plot(x=as.vector(ceiling(g_hic_matrix.norm.part.fix)),y=as.vector(ceiling(n_hic_matrix.raw.part)),pch=16,col="#5C88CC",xlim=c(0,200),ylim=c(0,200),xaxt="n",yaxt="n",xlab="",ylab="")
axis(side=1,at=seq(0,200,40),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=seq(0,200,40),y=-30 ,labels = seq(0,200,40),cex = 3,xpd=T,srt = 0)

axis(side=2,at=seq(0,200,40),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(y=seq(0,200,40),x=-30 ,labels = seq(0,200,40),cex = 3,xpd=T,srt = 0)
box(lwd=3)

cor(x=as.vector(ceiling(g_hic_matrix.norm.part.fix)),y=as.vector(ceiling(n_hic_matrix.raw.part)))

View(hela_ab_table[hela_ab_table$chrom_name=="chr14" & hela_ab_table$start >= 20e6 & hela_ab_table$start <= 30e6,])

chrom_NAIR_table = NAIR_table.merge[NAIR_table.merge$V1=="chr13",]


chrom_NAIR_table[which((chrom_NAIR_table$V3 - chrom_NAIR_table$V2)>1e6),]

# A/B compartment table 
binsize = 1e5
hela_ab_table = read.table("~/menghw_HD/data_table/Hela-genome-hic_AB_all_100000_ice.table",header = T,sep = "\t")

# A/B compartment length
HG19_LEN.compartment = nrow(hela_ab_table) * binsize 

# A/B ration in whole genome
a_ratio = sum(hela_ab_table$value[hela_ab_table$info=="A"]) / nrow(hela_ab_table)
b_ratio = sum(hela_ab_table$value[hela_ab_table$info=="B"]) / nrow(hela_ab_table)

##################################################################################
# fig4 ActD treatment: the enrichment changes suggests structure 
##################################################################################

##################################
# fig4 A all-all heatmap
##################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v02.RData")
load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-all-matrix_1000000.RData")
load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-n-hic-ActD_all_matrix_1000000.RData")
quantile(g_hic_all_matrix.norm,prob=0.9805)
quantile(n_hic_all_matrix,prob=0.985)
quantile(hela_n_hic_ActD_all_matrix,prob=0.968)

png(filename = "~/menghw_HD/R_image/Hela-hic-figure/Hela-g-hic_all_matrix.png",width = 4000,height = 4000)
par(mar=c(2,2,2,2))
plot.matrix(g_hic_all_matrix.norm,max_bound = 0.9805)
dev.off()

png(filename = "~/menghw_HD/R_image/Hela-hic-figure/Hela-n-hic_all_matrix.png",width = 4000,height = 4000)
par(mar=c(2,2,2,2))
plot.matrix(n_hic_all_matrix,max_bound = 0.985)
dev.off()

png(filename = "~/menghw_HD/R_image/Hela-hic-figure/Hela-n-hic-ActD_all_matrix.png",width = 4000,height = 4000)
par(mar=c(2,2,2,2))
plot.matrix(hela_n_hic_ActD_all_matrix,max_bound = 0.968)
dev.off()

###################################
# fig4 B chromosome changes
###################################
rm(list = ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v02.RData")
load("~/menghw_HD/data_table/fix_table/Hela-enrichment_v1.RData")


load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_track_chr_14_100000.RData")
n_hic_matrix_ActD.raw = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_100000/chr_14_100000_MAPQ20.txt",header = F,sep = ",")
n_hic_matrix_ActD.raw.part = matrix.part(n_hic_matrix_ActD.raw,chrom_index = 14,chrom_start = 19.1e6,binsize = 100e3)
g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = 14,chrom_start = 19.1e6,binsize = 100e3)
n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = 14,chrom_start = 19.1e6,binsize = 100e3)


load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_14_10000.RData")
n_hic_matrix_ActD.raw = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_10000/chr_14_10000_MAPQ20.txt",header = F,sep = ",")
n_hic_matrix_ActD.raw.part = matrix.part(n_hic_matrix_ActD.raw,chrom_index = 14,chrom_start = 20e6,chrom_end = 25e6,binsize = 10e3)
g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = 14,chrom_start = 20e6,chrom_end = 25e6,binsize = 10e3)
n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = 14,chrom_start = 20e6,chrom_end = 25e6,binsize = 10e3)


load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_13_10000.RData")
# n_hic_matrix_ActD.raw = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_10000/chr_13_10000_MAPQ20.txt",header = F,sep = ",")
# n_hic_matrix_ActD.raw.part = matrix.part(n_hic_matrix_ActD.raw,chrom_index = 14,chrom_start = 18e6,chrom_end = 21e6,binsize = 10e3)
g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = 13,chrom_start = 22e6,chrom_end = 25e6,binsize = 10e3)
n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = 13,chrom_start = 22e6,chrom_end = 25e6,binsize = 10e3)


par(mar=c(2,2,2,2))
max_bound_quantile = 0.95
quantile(g_hic_matrix.norm.part,prob= 0.95)
quantile(n_hic_matrix.raw.part,prob= 0.94)
quantile(n_hic_matrix_ActD.raw.part,prob= 0.992)

plot.matrix(g_hic_matrix.norm.part,max_bound = 0.96)
plot.matrix(n_hic_matrix.raw.part,max_bound = 0.92)
plot.matrix(n_hic_matrix_ActD.raw.part,max_bound = 0.992)


##################################################################################
# Fig 6 LAD NAD NAIR 
##################################################################################
load(file="~/menghw_HD/R_code/my_function/MyPlot_v02.RData")
png("~/menghw_HD/R_image/Hela-hic-figure/Hela-enrichment-region.png",width = 6000,height = 3000)


LAD_table = read.table("~/menghw_HD/data_table/fix_table/cLAD_human_fix.table",header = T,sep = "\t")

fig.facet <- layout(matrix(c(1:70),nrow = 70,byrow = FALSE),width = rep(10,70),heights = rep(1,70))
layout.show(fig.facet)

for(chrom_index in (1:23)){
  chrom_start = 0
  chrom_end = chrom.length(1)
  
  # plot NAD
  par(mar=c(0.5,1,0,1))
  plot.peak(df.peak = hela_n_seq_peak.merge,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#7080D7",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  # plot Hela-n-hic
  par(new=T)
  plot.peak(df.peak = hela_n_hic_peak.merge,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  # plot NAD NAIR overlap
  par(new=T)
  plot.peak(df.peak = hela_n_overlap,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#330674",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  # plot LAD
  plot.peak(df.peak = LAD_table,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "blue",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  # plot G-Band
  par(mar=c(0.8,1,0,1))
  plot.chromosome(chrom_index,chrom_start,chrom_end)
}
dev.off()


par(mar=c(3,3,3,3))
plot.chromosome(chrom_index=16,chrom_start = 0,chrom_end = chrom.length(chrom_index = 16),border.lwd = 5)
dev.off()

##################################################################################
# fig2 A Fixation
##################################################################################
rm(list=ls())
#####################################
# load data 
#####################################
load(file="~/menghw_HD/R_code/my_function/MyPlot_v05.RData")
# Hela_all_table 
hela_all_table = read.table("~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep = "\t")

# binsize table binsize=100K
hg19_loci_table = read.table("~/menghw_HD/data_table/hg19_loci_100000.bed",header = F,sep = "\t")
hg19_loci_table = cbind(hg19_loci_table,rep("genome_region",nrow(hg19_loci_table)))
colnames(hg19_loci_table) = c("chrom_name","start","end","name")

# NAD and NAIR
hela_n_seq_peak.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_w1e4.table",header = T,sep = "\t")
hela_n_hic_peak.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix_w1e5.table",header = F,sep = "\t")

hela_n_seq_peak.merge = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = T,sep = "\t")
hela_n_hic_peak.merge = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
colnames(hela_n_hic_peak.merge) = colnames(hela_n_seq_peak.merge)
hela_n_overlap = overlap(table_1 = hela_n_seq_peak.merge,table_2 = hela_n_hic_peak.merge,level_vetor = paste("chr",c(1:22,"X","Y"),sep = ""))

# RepeatMasker Data
load(file="~/menghw_HD/R_image/hg19_RepeatMasker_bed.RData")
SINE_coverage_table = read.table("~/menghw_HD/data_table/hg19_RM_SINE_coverage_100000.bed",header = F,sep = "\t")
colnames(SINE_coverage_table) = c("chrom_name","start","end","count","cover_len","binsize","rate")
LINE_coverage_table = read.table("~/menghw_HD/data_table/hg19_RM_LINE_coverage_100000.bed",header = F,sep = "\t")
colnames(LINE_coverage_table) = c("chrom_name","start","end","count","cover_len","binsize","rate")
LTR_coverage_table = read.table("~/menghw_HD/data_table/hg19_RM_LTR_coverage_100000.bed",header = F,sep = "\t")
colnames(LTR_coverage_table) = c("chrom_name","start","end","count","cover_len","binsize","rate")

LINE_LTR_coverage_table = cbind(LTR_coverage_table[,1:3],
                                LINE_coverage_table$count+LTR_coverage_table$count,
                                LINE_coverage_table$cover_len + LTR_coverage_table$cover_len,
                                LTR_coverage_table$binsize,
                                LTR_coverage_table$rate + LINE_coverage_table$rate)
colnames(LINE_LTR_coverage_table) = c("chrom_name","start","end","count","cover_len","binsize","rate")

# gene density table 
hg19_gene_table = read.table("~/menghw_HD/reference/gene_gtf/hg19_annotation.gene.bed",header = T,sep = "\t")
hg19_RNA_gene_cover_table = read.table("~/menghw_HD/reference/gene_gtf/hg19_RNA_gene_coverage_100000.bed",header = F,sep = "\t")
colnames(hg19_RNA_gene_cover_table) =  c("chrom_name","start","end","count","cover_len","binsize","rate")

hg19_no_RNA_gene_cover_table = read.table("~/menghw_HD/reference/gene_gtf/hg19_no_RNA_gene_coverage_100000.bed",header = F,sep = "\t")
colnames(hg19_no_RNA_gene_cover_table) = c("chrom_name","start","end","count","cover_len","binsize","rate")

# Histone Data coverage 
#### active
hela_H3K4me3_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_H3K4me3_ENCFF000BCO_coverage_100000.bed",header = F,sep = "\t")
hela_H3K36me3_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_H3K36me3_ENCFF000BCA_coverage_100000.bed",header = F,sep = "\t")
hela_POL2_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_POL2_ENCFF000PEL_coverage_100000.bed",header = F,sep = "\t")
#### repressive
hela_H3K27me3_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_H3K27me3_ENCFF000BBS_coverage_100000.bed",header = F,sep = "\t")
hela_H3K9me3_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_H3K9me3_ENCFF000BBG_coverage_100000.bed",header = F,sep = "\t")
hela_H4K20me3_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_H4K20me1_ENCFF000BDC_coverage_100000.bed",header = F,sep = "\t")

# coverage Data

# plot region
load("~/menghw_HD/R_code/my_function/MyPlot_v03.RData")
chrom_index = 17 
chrom_name = chrom.name(chrom_index)
chrom_start = 0
chrom_len = chrom.length(chrom_index)
chrom_end = chrom_len
chrom_binsize = 100000


#### make SINE & LINE
chrom_n_count = length(which(hela_all_table$chr_index==chrom_name & hela_all_table$g_cover_rate==0))

chrom_SINE_table = SINE_coverage_table[SINE_coverage_table$chrom_name==chrom_name,]
chrom_LINE_table = LINE_coverage_table[LINE_coverage_table$chrom_name==chrom_name,]
chrom_LINE_LTR_table = LINE_LTR_coverage_table[LINE_LTR_coverage_table$chrom_name==chrom_name,]

chrom_SINE_mean = sum(chrom_SINE_table$rate) / (nrow(chrom_SINE_table) - chrom_n_count)
chrom_LINE_mean = sum(chrom_LINE_table$rate) / (nrow(chrom_LINE_table) - chrom_n_count)
chrom_LINE_LTR_mean = sum(chrom_LINE_LTR_table$rate) / (nrow(chrom_LINE_LTR_table) - chrom_n_count)

chrom_SINE_table.fix = chrom_SINE_table[chrom_SINE_table$rate >= median(chrom_SINE_table$rate),]
chrom_LINE_table.fix = chrom_LINE_table[chrom_LINE_table$rate >= median(chrom_LINE_table$rate),]

chrom_SINE_table.fix = chrom_SINE_table
chrom_SINE_table.fix$rate = chrom_SINE_table$rate - chrom_SINE_mean
chrom_SINE_table.fix$rate = chrom_SINE_table$rate - mean(chrom_SINE_table$rate)
chrom_SINE_table.fix$rate[chrom_SINE_table$rate==0] = 0
chrom_SINE_table.fix = chrom_SINE_table.fix[chrom_SINE_table.fix$rate>=0,]

chrom_LINE_table.fix = chrom_LINE_table
chrom_LINE_table.fix$rate = chrom_LINE_table$rate - chrom_LINE_mean
chrom_LINE_table.fix$rate = chrom_LINE_table$rate - mean(chrom_LINE_table$rate)
chrom_LINE_table.fix$rate[chrom_LINE_table$rate==0] = 0
chrom_LINE_table.fix = chrom_LINE_table.fix[chrom_LINE_table.fix$rate>=0,]

chrom_LINE_LTR_table.fix = chrom_LINE_LTR_table
chrom_LINE_LTR_table.fix$rate = chrom_LINE_LTR_table$rate - chrom_LINE_LTR_mean
chrom_LINE_LTR_table.fix$rate = chrom_LINE_LTR_table$rate - mean(chrom_LINE_LTR_table$rate)
chrom_LINE_LTR_table.fix$rate[chrom_LINE_LTR_table$rate==0] = 0
chrom_LINE_LTR_table.fix = chrom_LINE_LTR_table.fix[chrom_LINE_LTR_table.fix$rate >=0,]

chrom_RNA_gene_table = hg19_RNA_gene_cover_table[hg19_RNA_gene_cover_table$chrom_name==chrom_name,]
chrom_RNA_gene_mean = sum(chrom_RNA_gene_table$rate) / (nrow(chrom_RNA_gene_table) - chrom_n_count)
chrom_RNA_gene_table.fix = chrom_RNA_gene_table[chrom_RNA_gene_table$rate>=chrom_RNA_gene_mean,]

chrom_no_RNA_gene_table = hg19_no_RNA_gene_cover_table[hg19_no_RNA_gene_cover_table$chrom_name==chrom_name,]
chrom_no_RNA_gene_mean = sum(chrom_no_RNA_gene_table$rate) / (nrow(chrom_no_RNA_gene_table) - chrom_n_count)
chrom_no_RNA_gene_table.fix = chrom_no_RNA_gene_table[chrom_no_RNA_gene_table$rate>=chrom_no_RNA_gene_mean,]
chrom_no_RNA_gene_table.fix = chrom_no_RNA_gene_table[chrom_no_RNA_gene_table$rate>=quantile(chrom_no_RNA_gene_table$rate,prob=0.8),]

chrom_H3K9me3_table = hela_H3K9me3_table[hela_H3K9me3_table$V1==chrom_name,]
chrom_H3K4me3_table = hela_H3K4me3_table[hela_H3K4me3_table$V1==chrom_name,]
chrom_H3K27me3_table = hela_H3K27me3_table[hela_H3K27me3_table$V1==chrom_name,]
chrom_H4K20me3_table = hela_H4K20me3_table[hela_H4K20me3_table$V1==chrom_name,]

# histone modification RPKM
chrom_H3K9me3_table$V7 = chrom_H3K9me3_table$V4 / sum(chrom_H3K9me3_table$V4) / (chrom_H3K9me3_table$V3 - chrom_H3K9me3_table$V2 + 1) * 1e9

chrom_H3K4me3_table$V7 = chrom_H3K4me3_table$V4 / sum(chrom_H3K4me3_table$V4) / (chrom_H3K4me3_table$V3 - chrom_H3K4me3_table$V2 + 1) * 1e9
chrom_H3K4me3_table.value = chrom_H3K4me3_table$V7 - median(chrom_H3K4me3_table$V7)
chrom_H3K4me3_table.value[chrom_H3K4me3_table$V7==0] = 0

chrom_H3K27me3_table$V7 = chrom_H3K27me3_table$V4 / sum(chrom_H3K27me3_table$V4) / (chrom_H3K27me3_table$V3 - chrom_H3K27me3_table$V2 + 1) * 1e9
chrom_H4K20me3_table$V7 = chrom_H4K20me3_table$V4 / sum(chrom_H4K20me3_table$V4) / (chrom_H4K20me3_table$V3 - chrom_H4K20me3_table$V2 + 1) * 1e9

fig.facet <- layout(matrix(c(1:7),nrow = 7,byrow = FALSE),width = rep(10,7),heights = rep(1,7))
layout.show(fig.facet)

#### cytoband plot
par(mar=c(3,5,3,1))
plot.chromosome(chrom_index,chrom_start,chrom_end)
dev.off()
#### SINE & LINE plot 
# plot(x=chrom_SINE_table$start,y=chrom_SINE_table$rate,type="h")
# plot(x=chrom_LINE_table$start,y=chrom_LINE_table$rate,type="h")

# par(mar=c(0.1,1,0,1))
par(mar=c(3,5,3,1),family="Arial",font=2)
plot(x=chrom_SINE_table.fix$start,y=chrom_SINE_table.fix$rate,type="h",ylim=c(0,0.3),lwd=1,xlim=c(chrom_start,chrom_end),xaxt="n",yaxt="n",frame.plot=F,col="#012E34",xlab="",ylab="")
axis(side=2,at=seq(0,0.3,0.1),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=par()$usr[1] - 300000,y=seq(0,0.3,0.1),labels = c(0,0.1,0.2,0.3),cex = 3.5,xpd=T,pos=2)


# par(mar=c(0.5,1,0,1))
# plot(x=chrom_LINE_table.fix$start,y=chrom_LINE_table.fix$rate,type="h",ylim=c(-0.3,0.5),lwd=1,xlim=c(chrom_start,chrom_end),xaxt="n",yaxt="n",frame.plot=F)

par(mar=c(0.25,1,0,1))
plot(x=chrom_LINE_LTR_table.fix$start,y=chrom_LINE_LTR_table.fix$rate,type="h",ylim=c(0,0.3),lwd=1,xlim=c(chrom_start,chrom_end),xaxt="n",yaxt="n",frame.plot=F,col="#805E15",ylab="")
axis(side=2,at=seq(0,0.3,0.1),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=par()$usr[1] - 300000,y=seq(0,0.3,0.1),labels = c(0,0.1,0.2,0.3),cex = 3.5,xpd=T,pos=2)

#### gene density
par(mar=c(0.25,1,0,1))
plot(x=chrom_RNA_gene_table.fix$start,y=chrom_RNA_gene_table.fix$rate,type="h",ylim=c(0,1),lwd=0.5,xlim=c(chrom_start,chrom_end),xaxt="n",yaxt="n",frame.plot=F,col="#0E464E",ylab="")
axis(side=2,at=seq(0,1,0.5),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=par()$usr[1] - 300000,y=seq(0,1,0.5),labels = c(0,0.5,1),cex = 3.5,xpd=T,pos=2)

par(mar=c(0.25,1,0,1))
plot(x=chrom_no_RNA_gene_table.fix$start,y=chrom_no_RNA_gene_table.fix$rate,type="h",ylim=c(0,1),lwd=0.5,xlim=c(chrom_start,chrom_end),xaxt="n",yaxt="n",frame.plot=F,col="#69969C")


#### NAIR plot
par(mar=c(0.25,1,0,1))
plot.peak(df.peak = hela_n_seq_peak.merge,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#7080D7",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

# plot Hela-n-hic
par(new=T)
plot.peak(df.peak = hela_n_hic_peak.merge,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

# plot NAD NAIR overlap
par(new=T)
plot.peak(df.peak = hela_n_overlap,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#330674",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)


#### histone modification 
par(mar=c(0.25,1,0,1))
plot(x=chrom_H3K4me3_table$V2,y=chrom_H3K4me3_table.value,type="h",ylim=c(-5,10),lwd=0.5,xlim=c(chrom_start,chrom_end),xaxt="n",yaxt="n",frame.plot=F,ylab="")
axis(side=2,at=seq(-5,10,5),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=par()$usr[1] - 300000,y=seq(-5,10,5),labels = seq(-5,10,5),cex = 3.5,xpd=T,pos=2)

# plot(x=chrom_H3K9me3_table$V2,y=chrom_H3K9me3_table$V7-median(chrom_H3K9me3_table$V7),type="h",ylim=c(0,10),lwd=0.5,xlim=c(chrom_start,chrom_end),xaxt="n",yaxt="n",frame.plot=F)

# plot(x=chrom_H3K27me3_table$V2,y=chrom_H3K27me3_table$V7-median(chrom_H3K27me3_table$V7),type="h",ylim=c(0,10),lwd=0.5,xlim=c(chrom_start,chrom_end),xaxt="n",yaxt="n",frame.plot=F)

# plot(x=chrom_H4K20me3_table$V2,y=chrom_H4K20me3_table$V7-median(chrom_H4K20me3_table$V7),type="h",ylim=c(0,10),lwd=0.5,xlim=c(chrom_start,chrom_end),xaxt="n",yaxt="n",frame.plot=F)

chrom_hela_g_hic_TAD_table = read.table("~/menghw_HD/Project/TAD_calling/TAD_region/Hela-g-hic_chr_13_20000_65.tad",header = F,sep = "\t")
chrom_hela_g_hic_TAD_table = cbind(chrom_hela_g_hic_TAD_table,rep(1,nrow(chrom_hela_g_hic_TAD_table)),rep(1,nrow(chrom_hela_g_hic_TAD_table)))
colnames(chrom_hela_g_hic_TAD_table) = c("chrom_name","start","end","value-1","value-2")

chrom_hela_n_hic_TAD_table = read.table("~/menghw_HD/Project/TAD_calling/TAD_region/Hela-n-hic_chr_13_20000_65.tad",header = F,sep = "\t")
chrom_hela_n_hic_TAD_table = cbind(chrom_hela_n_hic_TAD_table,rep(1,nrow(chrom_hela_n_hic_TAD_table)))
colnames(chrom_hela_n_hic_TAD_table) = c("chrom_name","start","end","value")

load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_13_20000.RData")
g_hic_matrix.norm.part = matrix.part(g_hic_matrix.norm,chrom_index = 13,chrom_start = 20e6,chrom_end = 30e6,binsize = 20e3)
n_hic_matrix.raw.part = matrix.part(n_hic_matrix,chrom_index = 13,chrom_start = 20e6,chrom_end = 30e6,binsize = 20e3)


par(mar=c(2,2,2,2))
max_bound_quantile = 0.95
quantile(g_hic_matrix.norm.part,prob= 0.97)
quantile(g_hic_matrix.norm.part.fix,prob= 0.95)
quantile(n_hic_matrix.raw.part,prob= 0.95)
quantile(n_hic_matrix_ActD.raw.part,prob= 0.96)


quantile(log2(g_hic_matrix.norm.part+1),prob= 0.976)
quantile(log2(n_hic_matrix.raw.part+1),prob= 0.98)

##################################################################################
# NAD-NAD NAT-NAT NAT-other 比例
##################################################################################
chrom.index <- function(chrom_name = "chr1"){
  chrom_index = 0
  if(chrom_name == "rDNA"){
    chrom_index = 25
  }else if(chrom_name == "chrY"){
    chrom_index = 24
  }else if(chrom_name == "chrX"){
    chrom_index = 23
  }else{
    chrom_index = as.numeric(substr(chrom_name,4,nchar(chrom_name)))
  }
  return(chrom_index)
}

get.genome_index <- function(chrom_index,chrom_postion=0,binsize=1e6){
  acc_index = 0
  if(chrom_index >= 2){
    for(i in c(1:(chrom_index-1))){
      acc_index.part = (chrom.length(i) %/% binsize) + 1
      acc_index = acc_index + acc_index.part
    }
  }
  chrom_bin_index = (chrom_postion %/% binsize) + 1
  return(acc_index + chrom_bin_index)
}

# load function and data
load(file="~/menghw_HD/R_code/my_function/MyPlot_v04.RData")
load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-all-matrix_1000000.RData")

hela_NAT = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
colnames(hela_NAT) = c("chrom_name","start","end","name","value")
hela_NAD = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep = "\t")
colnames(hela_NAD) = c("chrom_name","start","end","name","value")

binsize = 1e6
####################################
## NAT add genome index 
####################################
NAT_start_index_list = NULL
NAT_end_index_list = NULL
for(i in c(1:nrow(hela_NAT))){
  NAT_line = hela_NAT[i,]
  NAT_chrom_index = chrom.index(chrom_name = as.character(NAT_line$chrom_name))
  NAT_start_index = get.genome_index(chrom_index = NAT_chrom_index,chrom_postion = NAT_line$start,binsize)
  NAT_end_index = get.genome_index(chrom_index = NAT_chrom_index,chrom_postion = NAT_line$end,binsize)
  NAT_start_index_list = c(NAT_start_index_list,NAT_start_index)
  NAT_end_index_list = c(NAT_end_index_list,NAT_end_index)
}


hela_NAT.fix = cbind(hela_NAT,NAT_start_index_list,NAT_end_index_list)
colnames(hela_NAT.fix) = c("chrom_name","start","end","name","value","start_index","end_index")
# write.table(hela_NAT.fix,file="~/menghw_HD/data_table/Hela_NAT_genome_index.bed",col.names = F,row.names = F,quote = F,sep = "\t")

####################################
## NAD add genome index 
####################################
NAD_start_index_list = NULL
NAD_end_index_list = NULL
for(i in c(1:nrow(hela_NAD))){
  NAD_line = hela_NAD[i,]
  NAD_chrom_index = chrom.index(chrom_name = as.character(NAD_line$chrom_name))
  NAD_start_index = get.genome_index(chrom_index = NAD_chrom_index,chrom_postion = NAD_line$start,binsize)
  NAD_end_index = get.genome_index(chrom_index = NAD_chrom_index,chrom_postion = NAD_line$end,binsize)
  NAD_start_index_list = c(NAD_start_index_list,NAD_start_index)
  NAD_end_index_list = c(NAD_end_index_list,NAD_end_index)
}

hela_NAD.fix = cbind(hela_NAD,NAD_start_index_list,NAD_end_index_list)
colnames(hela_NAD.fix) = c("chrom_name","start","end","name","value","start_index","end_index")

####################################
## make NAD-NAD matrix 
####################################
# load genome all-all matrix 
all_matrix = g_hic_all_matrix
all_matrix = n_hic_all_matrix

NAD_matrix.raw = matrix(0,nrow =nrow(hela_NAD.fix),ncol = nrow(hela_NAD.fix) )

for(i in c(1:nrow(hela_NAD.fix))){
  row.start_index = hela_NAD.fix$start_index[i]
  row.end_index = hela_NAD.fix$end_index[i]
  for(j in c(i:nrow(hela_NAD.fix))){
    if(i == j){
      inter_cis.matrix = all_matrix[row.start_index:row.end_index,row.start_index:row.end_index]
      inter_cis =sum(inter_cis.matrix[upper.tri(inter_cis.matrix,diag = T)])
      NAD_matrix.raw[i,j] = inter_cis
    }else{
      col.start_index = hela_NAD.fix$start_index[j]
      col.end_index = hela_NAD.fix$end_index[j]
      inter_trans.matrix = all_matrix[row.start_index:row.end_index,col.start_index:col.end_index]
      inter_trans = sum(inter_trans.matrix)
      NAD_matrix.raw[i,j] = inter_trans
      NAD_matrix.raw[j,i] = inter_trans
    }
  }
}

total_inter = sum(all_matrix[upper.tri(all_matrix,diag = T)])
NAD_inter = sum(NAD_matrix.raw[upper.tri(NAD_matrix.raw,diag = T)])

g_hic_info_NAD = c(total_inter,NAD_inter,total_inter - NAD_inter)
n_hic_info_NAD = c(total_inter,NAD_inter,total_inter - NAD_inter)


NAD_matrix.fix = NAD_matrix.raw / rep((hela_NAD.fix$end - hela_NAD.fix$start),each=nrow(hela_NAD.fix))* 1e6
plot.matrix(log2(NAD_matrix.fix+1),max_bound = 1)
plot.matrix(log2(NAD_matrix.raw+1),max_bound = 1)

####################################
## make NAT-NAT matrix 
####################################
NAT_matrix.raw = matrix(0,nrow =nrow(hela_NAT.fix),ncol = nrow(hela_NAT.fix) )

for(i in c(1:nrow(hela_NAT.fix))){
  row.start_index = hela_NAT.fix$start_index[i]
  row.end_index = hela_NAT.fix$end_index[i]
  for(j in c(i:nrow(hela_NAT.fix))){
    if(i == j){
      inter_cis.matrix = all_matrix[row.start_index:row.end_index,row.start_index:row.end_index]
      inter_cis =sum(inter_cis.matrix[upper.tri(inter_cis.matrix,diag = T)])
      NAT_matrix.raw[i,j] = inter_cis
    }else{
      col.start_index = hela_NAT.fix$start_index[j]
      col.end_index = hela_NAT.fix$end_index[j]
      inter_trans.matrix = all_matrix[row.start_index:row.end_index,col.start_index:col.end_index]
      inter_trans = sum(inter_trans.matrix)
      NAT_matrix.raw[i,j] = inter_trans
      NAT_matrix.raw[j,i] = inter_trans
    }
  }
}

# total_inter = sum(all_matrix[upper.tri(all_matrix,diag = T)])
# NAT_inter = sum(NAT_matrix.raw[upper.tri(NAT_matrix.raw,diag = T)])

# g_hic_info_NAT = c(total_inter,NAT_inter,total_inter-NAT_inter)
# n_hic_info_NAT = c(total_inter,NAT_inter,total_inter-NAT_inter)



NAT_matrix.fix = NAT_matrix.raw / rep((hela_NAT.fix$end - hela_NAT.fix$start),each=nrow(NAT_matrix.raw)) * 1e6
# hela_g_hic_NAT_matrix.raw = NAT_matrix.raw
# hela_g_hic_NAT_matrix.fix = NAT_matrix.raw / rep((hela_NAT.fix$end - hela_NAT.fix$start),each=nrow(NAT_matrix.raw)) * 1e6

hela_n_hic_NAT_matrix.raw = NAT_matrix.raw
hela_n_hic_NAT_matrix.fix = NAT_matrix.raw / rep((hela_NAT.fix$end - hela_NAT.fix$start),each=nrow(NAT_matrix.raw)) * 1e6

plot.matrix(log2(NAT_matrix.fix+1),max_bound = 1,min_bound=0.5)
plot.matrix((NAT_matrix.fix),max_bound = 0.95)

plot.matrix(log2(hela_n_hic_NAT_matrix.fix+1),max_bound = 1,min_bound=0)
plot.matrix((hela_n_hic_NAT_matrix.fix),max_bound = 0.95)

root_matrix = sqrt(colSums(hela_n_hic_NAT_matrix.raw)) %*% t(sqrt(colSums(hela_n_hic_NAT_matrix.raw)))
hela_n_hic_NAT_matrix.fix.fix =  hela_n_hic_NAT_matrix.raw / root_matrix

root_matrix = sqrt(colSums(hela_g_hic_NAT_matrix.raw)) %*% t(sqrt(colSums(hela_g_hic_NAT_matrix.raw)))
hela_g_hic_NAT_matrix.fix.fix =  hela_g_hic_NAT_matrix.raw / root_matrix

plot.matrix(hela_g_hic_NAT_matrix.fix.fix,max_bound = 0.975)
plot.matrix(hela_n_hic_NAT_matrix.fix.fix,max_bound = 0.965)


# barplot.mat => 0.3447659 0.6552341 0.6072660 0.3927340
barplot.mat = matrix(1,2,2)
barplot.mat[1,1] = g_hic_info_NAT[2] / g_hic_info_NAT[1]
barplot.mat[2,1] = g_hic_info_NAT[3] / g_hic_info_NAT[1]

barplot.mat[1,2] = n_hic_info_NAT[2] / n_hic_info_NAT[1]
barplot.mat[2,2] = n_hic_info_NAT[3] / n_hic_info_NAT[1]

####################################
## NAD-NAD interaction barplot 
####################################
par(mar=c(11,8,2,1),family="Arial",font=1)
barplot(barplot.mat,axes = F,col=c("#FF7304","#5CCDC9"),width = 0.5,space = 0.25)
axis(side=2,at=seq(0,1,0.25),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=-0.1,y=seq(0,1,0.25),labels = paste0(seq(0,100,25),"%"),cex = 3,xpd=T,pos = NULL)
text(x=c(0.45,1.05),y=-0.05,labels = c("In Situ Hi-C","N-Hi-C"),cex = 3,xpd=T,srt = 45,pos = 2)


########################################################################
## data correlation
########################################################################
# hela_g_seq, hela_g_hic,hela_n_seq,hela_n_hic
relation_matrix = matrix(rep(1,16),4)
relation_matrix[1,2] = cor(x=hela_all_table$g_count,y=hela_all_table$g_hic)
relation_matrix[2,1] = relation_matrix[1,2]

relation_matrix[1,3] = cor(x=hela_all_table$g_count,y=hela_all_table$n_count)
relation_matrix[3,1] = relation_matrix[1,3]

relation_matrix[1,4] = cor(x=hela_all_table$g_count,y=hela_all_table$n_hic)
relation_matrix[4,1] = relation_matrix[1,4]

relation_matrix[2,3] = cor(x=hela_all_table$g_hic,y=hela_all_table$n_count)
relation_matrix[3,2] = relation_matrix[2,3]


relation_matrix[2,4] = cor(x=hela_all_table$g_hic,y=hela_all_table$n_hic)
relation_matrix[4,2] = relation_matrix[2,3]

relation_matrix[3,4] = cor(x=hela_all_table$n_count,y=hela_all_table$n_hic)
relation_matrix[4,3] = relation_matrix[3,4]

par(mar=c(1,1,1,1))
plot.matrix(relation_matrix-0.1,max_bound = 1)


##################################################################################
# fig4 heatmap-track-heatmap
##################################################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")

###############################
## load matrix data 
###############################
chrom_index = 4

load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_4_100000.RData")
hela_ab_value_table = read.table("~/menghw_HD/data_table/Hela_AB_compart_value.table",header = T,sep = "\t")
hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
hela_NAD_table = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep = "\t")

chrom_name = chrom.name(chrom_index)
chrom_start = 0
chrom_end = chrom.length(chrom_index)

fig.facet <- layout(matrix(c(1:5),nrow = 5,ncol = 1),width = rep(10,5),heights = c(10,1.5,1,1,10))
layout.show(fig.facet)

###############################
## first column
###############################

# hela-g-hic matrix
par(mar=c(0.5,2,0.5,2),family="Arial",font=1)
plot.matrix(g_hic_matrix.norm,bound.max = 0.95)

# hela-g-hic.norm A/B compartment
par(mar=c(0.5,2,0,2),family="Arial",font=1)
ab_x = start(pc.norm$PC1)
ab_y = score(pc.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",yaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))

# Hela-NAD
par(mar=c(0,2,0,2),family="Arial",font=1)
plot.peak(df.peak = hela_NAD_table,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#7080D7",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)


# Hela-NAIR
par(mar=c(0,2,0,2),family="Arial",font=1)
plot.peak(df.peak = hela_NAT_table,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)


# hela-n-hic matrix
par(mar=c(0.5,2,0,2),family="Arial",font=1)
plot.matrix(n_hic_matrix,bound.max = 0.95)


plot.chromosome(chrom_index = 4,chrom_start = 0,chrom_end = chrom_end)

sum(hela_NAD_table$V3 - hela_NAD_table$V2)
sum(hela_NAT_table$V3 - hela_NAT_table$V2)


dev.off()
























