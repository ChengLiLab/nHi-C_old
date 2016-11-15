################################################################
# fig6 rDNA 富集倍数 与cis trans 改变比例
################################################################
rm(list=ls())
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")

hic_count.hela_g_hic.rDNA.cis = 200979
hic_count.hela_g_hic.rDNA.trans = 162873

hic_count.hela_n_hic.rDNA.cis = 2478769
hic_count.hela_n_hic.rDNA.trans = 637220

hic_count.hela_g_hic_ActD.rDNA.cis = 69103
hic_count.hela_g_hic_ActD.rDNA.trans = 181594

hic_count.hela_n_hic_ActD.rDNA.cis = 491898
hic_count.hela_n_hic_ActD.rDNA.trans = 993802

hic_count_table = data.frame(hela_g_hic = c(hic_count.hela_g_hic.rDNA.cis,hic_count.hela_g_hic.rDNA.trans),
                             hela_g_hic_ActD = c(hic_count.hela_g_hic_ActD.rDNA.cis,hic_count.hela_g_hic_ActD.rDNA.trans),
                             hela_n_hic = c(hic_count.hela_n_hic.rDNA.cis,hic_count.hela_n_hic.rDNA.trans),
                             hela_n_hic_ActD = c(hic_count.hela_n_hic_ActD.rDNA.cis,hic_count.hela_n_hic_ActD.rDNA.trans))


hic_count_table = data.frame(hela_g_hic = c(hic_count.hela_g_hic.cis,hic_count.hela_g_hic.trans),
                             hela_g_hic_ActD = c(hic_count.hela_g_hic_ActD.cis,hic_count.hela_g_hic_ActD.trans),
                             hela_n_hic = c(hic_count.hela_n_hic.cis,hic_count.hela_n_hic.trans),
                             hela_n_hic_ActD = c(hic_count.hela_n_hic_ActD.cis,hic_count.hela_n_hic_ActD.trans))
rownames(hic_count_table) = c("Cis_count","Trans_count")

barplot.matrix = rbind(as.matrix(hic_count_table)[1,] / colSums(as.matrix(hic_count_table)),as.matrix(hic_count_table)[2,] / colSums(as.matrix(hic_count_table)))
barplot.matrix = as.matrix(hic_count_table) / c(hic_count.hela_g_hic.total,hic_count.hela_g_hic_ActD.total,hic_count.hela_n_hic.total,hic_count.hela_n_hic_ActD.total) * 300e6

# par(mar=c(8,6,2,1),family="Arial",font=1)
par(mar=c(2,2,1,1),family="Arial",font=1)
barplot(barplot.matrix,col = c("#C936D3","#89A110"),border = F,axes = F,axisnames=F)
axis(side=2,at=seq(0,1,0.25),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=-0.9,y=seq(0,1,0.25),labels = paste0(seq(0,100,25),"%"),cex = 2,xpd=T,pos = NULL)
text(x=c(1.2,2.2,3.4,4.6),y=-0.05,labels = c("Hi-C","Hi-C ActD","nHi-C","nHi-C ActD"),cex = 2,xpd=T,srt = 45,pos = 2)

# plot rDNA total interaction 
hic_count_table = data.frame(hela_g_hic = c(hic_count.hela_g_hic.rDNA.cis,hic_count.hela_g_hic.rDNA.trans),
                             hela_g_hic_ActD = c(hic_count.hela_g_hic_ActD.rDNA.cis,hic_count.hela_g_hic_ActD.rDNA.trans),
                             hela_n_hic = c(hic_count.hela_n_hic.rDNA.cis,hic_count.hela_n_hic.rDNA.trans),
                             hela_n_hic_ActD = c(hic_count.hela_n_hic_ActD.rDNA.cis,hic_count.hela_n_hic_ActD.rDNA.trans))
rownames(hic_count_table) = c("Cis_count","Trans_count")

barplot.matrix = colSums(as.matrix(hic_count_table)) / c(hic_count.hela_g_hic.total,hic_count.hela_g_hic_ActD.total,hic_count.hela_n_hic.total,hic_count.hela_n_hic_ActD.total) * 300e2
par(mar=c(2,2,1,1),family="Arial",font=1)
barplot(barplot.matrix,col = c("#D4BA6A","#D4BA6A","#565695","#565695"),border = F,axes = F,axisnames=F,ylim=c(0,300))
axis(side=2,at=seq(0,300,100),labels = F,cex.axis=3,lwd=3,tck=-0.05)


################################################################
# fig6 plot heatmap
################################################################
rm(list=ls())

hela_n_hic_rDNA.matrix = read.table("~/menghw_HD/our_data/rDNA_data/Hela-n-hic_chr_25_1000.txt",header = F,sep = ",")
hela_n_hic_ActD_rDNA.matrix = read.table("~/menghw_HD/our_data/rDNA_data/Hela-n-hic-ActD_chr_25_1000.txt",header = F,sep = ",")
load(file = "~/menghw_HD/R_code/my_function/MyPlot_v08.RData")


par(mar=c(1,1,1,1))
plot.matrix(hela_n_hic_rDNA.matrix,bound.max = 0.95)
plot.matrix(hela_n_hic_ActD_rDNA.matrix,bound.max = 0.99)

##### fix the rDNA matrix with sequencing depth 
hela_n_hic_rDNA.matrix.fix = as.matrix(hela_n_hic_rDNA.matrix / (hic_count.hela_n_hic.total / hic_count.hela_n_hic_ActD.total))
hela_n_hic_ActD_rDNA.matrix.fix = as.matrix(hela_n_hic_ActD_rDNA.matrix)

hela_n_hic_rDNA.matrix.fix = round(as.matrix(hela_n_hic_rDNA.matrix / hic_count.hela_n_hic.total) * 300e6)
hela_n_hic_ActD_rDNA.matrix.fix = round(as.matrix(hela_n_hic_ActD_rDNA.matrix) / (hic_count.hela_g_hic_ActD.total) * 300e6)

hela_n_hic_rDNA.matrix.fix.fixmax = hela_n_hic_rDNA.matrix.fix
hela_n_hic_rDNA.matrix.fix.fixmax[hela_n_hic_rDNA.matrix.fix.fixmax>=quantile(hela_n_hic_rDNA.matrix.fix,0.86)] = round(quantile(hela_n_hic_rDNA.matrix.fix,0.86))
hela_n_hic_ActD_rDNA.matrix.fix.fixmax = hela_n_hic_ActD_rDNA.matrix.fix
hela_n_hic_ActD_rDNA.matrix.fix.fixmax[hela_n_hic_ActD_rDNA.matrix.fix.fixmax>= quantile(hela_n_hic_ActD_rDNA.matrix.fix,0.95)] = round(quantile(hela_n_hic_ActD_rDNA.matrix.fix,0.95))

# write table 
write.table(hela_n_hic_rDNA.matrix.fix.fixmax,file="~/menghw_HD/our_data/rDNA_data/Hela-n-hic_chr_25_1000_fixmax.txt",col.names = F,row.names = F,sep = "\t",quote = F)
write.table(hela_n_hic_ActD_rDNA.matrix.fix.fixmax,file="~/menghw_HD/our_data/rDNA_data/Hela-n-hic-ActD_chr_25_1000_fixmax.txt",col.names = F,row.names = F,sep = "\t",quote = F)


hela_n_hic_rDNA_differ.matrix =  hela_n_hic_ActD_rDNA.matrix.fix / hela_n_hic_rDNA.matrix.fix
hela_n_hic_rDNA_differ.matrix.log = log2(hela_n_hic_rDNA_differ.matrix)
hela_n_hic_rDNA_differ.matrix.log[hela_n_hic_ActD_rDNA.matrix.fix==0 | hela_n_hic_rDNA.matrix.fix==0] = 0

quantile(hela_n_hic_rDNA.matrix.fix,prob=seq(0.85,1,0.01))
quantile(hela_n_hic_ActD_rDNA.matrix.fix,prob=seq(0.95,1,0.01))
quantile(hela_n_hic_rDNA_differ.matrix.log,prob=c(0.05,0.95))

plot.matrix(hela_n_hic_rDNA.matrix.fix,bound.max = 0.86,n_block_color = "gray")
plot.matrix(hela_n_hic_ActD_rDNA.matrix.fix,bound.max = 0.95,n_block_color = "gray")
plot.matrix(hela_n_hic_rDNA_differ.matrix.log,bound.min = 0.00,col.min = "blue",col.max = "red",col.boundary = 0,bound.max = 0.9,n_block_color = "gray")


################################################################
# fig6 active enhancer
################################################################
filter = 1274440 / 2 * 100 / 2.8e9
hela_enhancer = read.table("~/menghw_HD/data_table/Hela_active_enhancer/Hela_active_enhancer_coverage.bed",header = F,sep = "\t")
hela_enhancer$V1 = factor(x=hela_enhancer$V1,levels = c(paste0("chr",c(1:22)),"chrX"))

hela_enhancer.active = hela_enhancer[hela_enhancer$V13 > filter * 10 & hela_enhancer$V12 > 1000,]

hela_enhancer.active.fix = hela_enhancer.active[,c(1,2,3,11,12,13)]
colnames(hela_enhancer.active.fix) = c("chrom_name","start","end","cover_length","cover_rate","length")

write.table(hela_enhancer.active.fix,file = "~/menghw_HD/data_table/Hela_active_enhancer/Hela_rDNA_enhancer.bed",col.names = T,row.names = F,sep = "\t",quote = F)


################################################################
# fig6 rDNA interaction with genome
################################################################
rm(list=ls())

load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")

hela_TR_count = read.table("~/menghw_HD/our_data/rDNA_data/rDNA_coverage/Hela-n-hic_rDNA_coverage_w1e5_TR.bed",header = F,sep = "\t")
hela_IGS_count = read.table("~/menghw_HD/our_data/rDNA_data/rDNA_coverage/Hela-n-hic_rDNA_coverage_w1e5_IGS.bed",header = F,sep = "\t")
hela_ActD_TR_count = read.table("~/menghw_HD/our_data/rDNA_data/rDNA_coverage/Hela-n-hic-ActD_rDNA_coverage_w1e5_TR.bed",header = F,sep = "\t")
hela_ActD_IGS_count = read.table("~/menghw_HD/our_data/rDNA_data/rDNA_coverage/Hela-n-hic-ActD_rDNA_coverage_w1e5_IGS.bed",header = F,sep = "\t")

colnames(hela_TR_count) = c("chrom_name","start","end","count","coverage_length","binsize","coverage_rate")
colnames(hela_IGS_count) = c("chrom_name","start","end","count","coverage_length","binsize","coverage_rate")
colnames(hela_ActD_TR_count) = c("chrom_name","start","end","count","coverage_length","binsize","coverage_rate")
colnames(hela_ActD_IGS_count) = c("chrom_name","start","end","count","coverage_length","binsize","coverage_rate")

sum(hela_TR_count$count) / sum(hela_IGS_count$count)
sum(hela_ActD_TR_count$count) / sum(hela_ActD_IGS_count$count)


quantile(hela_TR_count.vector,prob=seq(0.95,1,0.001))
plot(density(hela_TR_count.vector))
length(hela_TR_count.vector)  

# View(hela_TR_count[hela_TR_count$count>30 & hela_TR_count$count<=1000,])
table(hela_TR_count[hela_TR_count$count>0 & hela_TR_count$count<=1000,1])
table(hela_IGS_count[hela_IGS_count$count>0 & hela_IGS_count$count<=1000,1])

quantile(hela_IGS_count.vector,prob=seq(0.75,1,0.01))
quantile(hela_ActD_TR_count.vector,prob=seq(0.75,1,0.01))
quantile(hela_ActD_IGS_count.vector,prob=seq(0.75,1,0.01))

hela_TR_count.vector = hela_TR_count$count[hela_TR_count$count>1 & hela_TR_count$count<=1000]
hela_IGS_count.vector = hela_IGS_count$count[hela_IGS_count$count>1 & hela_IGS_count$count <= 1000]
hela_ActD_IGS_count.vector = hela_ActD_IGS_count$count[hela_ActD_IGS_count$count>1 & hela_ActD_IGS_count$count <= 1000]
hela_ActD_TR_count.vector = hela_ActD_TR_count$count[hela_ActD_TR_count$count>1 & hela_ActD_TR_count$count <=1000]

wilcox.test(hela_TR_count.vector,hela_IGS_count.vector)
wilcox.test(hela_ActD_IGS_count.vector,hela_ActD_TR_count.vector)

length(which((hela_ActD_IGS_count.vector > 100)))
length(which(hela_IGS_count.vector>100))

boxplot(log2(hela_TR_count.vector),log2(hela_IGS_count.vector),log2(hela_ActD_TR_count.vector),log2(hela_ActD_IGS_count.vector),ylim=c(0,10),col = c("#D4BA6A","#D4BA6A","#565695","#565695"),xaxt="n",yaxt="n",frame.plot=F)
axis(side=1,at=c(1:4),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
axis(side=2,at=seq(0,10,2),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)

################################################################
# fig6 rDNA interaction distance 
################################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/make_dist_table.RData")
load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")
rDNA_matrix = as.matrix(read.table("~/menghw_HD/our_data/rDNA_data/Hela-n-hic_chr_25_1000.txt",header = F,sep = ","))
rDNA_matrix.ActD = as.matrix(read.table("~/menghw_HD/our_data/rDNA_data/Hela-n-hic-ActD_chr_25_1000.txt",header = F,sep = ","))

rDNA_matrix.ActD.fix = round(rDNA_matrix.ActD / sum(rDNA_matrix.ActD) * abs(sum(rDNA_matrix) - sum(rDNA_matrix.ActD)))

chrom_index = 25
rDNA_dist_table = make.dist.table(rDNA_matrix,1000)
rDNA_ActD_dist_table = make.dist.table(rDNA_matrix.ActD,1000)
rDNA_ActD_dist_table = make.dist.table(rDNA_matrix.ActD.fix,1000)

rDNA_differ = rDNA_matrix.ActD.fix / rDNA_matrix
rDNA_differ[rDNA_matrix==0] = 0
rDNA_differ.log = log2(rDNA_differ)
rDNA_differ.log[rDNA_matrix==0 | rDNA_differ==0] = 0


par(mar=c(1,1,1,1))
plot.matrix(rDNA_matrix,bound.max = 0.93,n_block_color = "gray")
plot.matrix(rDNA_matrix.ActD.fix,bound.max = 0.945,n_block_color = "gray")
plot.matrix(rDNA_differ.log,bound.max = 0.95,col.min = "blue",bound.min = 0.05,col.boundary = 0,n_block_color = "gray")


quantile(rDNA_matrix,0.93)
quantile(rDNA_matrix.ActD.fix,0.947)
quantile(rDNA_differ.log,prob=c(0.05,0.95))

ylim = c(2,6)
xlim = c(0,chrom.length(chrom_index))

plot(x=(rDNA_dist_table$start),y=log10(rDNA_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#1533AD",frame.plot=F,yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
par(new=T)
plot(x=(rDNA_ActD_dist_table$start),y=log10(rDNA_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#AA3355",frame.plot=F,yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
axis(side=2,at=c(3:6))


plot(x=log10(rDNA_dist_table$start),y=log10(rDNA_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#1533AD",frame.plot=F,yaxt="n",xlab="",ylab="",ylim=ylim)
par(new=T)
plot(x=log10(rDNA_ActD_dist_table$start),y=log10(rDNA_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#AA3355",frame.plot=F,yaxt="n",xlab="",ylab="",ylim=ylim)
axis(side=2,at=c(3:6))




################################################################
# fig6 rDNA hot spot
################################################################
rm(list=ls())

load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")

# hela_TR_count = read.table("~/menghw_HD/our_data/rDNA_data/rDNA_coverage/Hela-n-hic_rDNA_coverage_w1e5_TR.bed",header = F,sep = "\t")
# hela_IGS_count = read.table("~/menghw_HD/our_data/rDNA_data/rDNA_coverage/Hela-n-hic_rDNA_coverage_w1e5_IGS.bed",header = F,sep = "\t")
# hela_ActD_TR_count = read.table("~/menghw_HD/our_data/rDNA_data/rDNA_coverage/Hela-n-hic-ActD_rDNA_coverage_w1e5_TR.bed",header = F,sep = "\t")
# hela_ActD_IGS_count = read.table("~/menghw_HD/our_data/rDNA_data/rDNA_coverage/Hela-n-hic-ActD_rDNA_coverage_w1e5_IGS.bed",header = F,sep = "\t")

hela_rDNA_trans_table = read.table("~/menghw_HD/our_data/rDNA_data/rDNA_coverage/Hela-n-hic_w1e5_coverage.bed",header = F,sep = "\t")
hela_rDNA_ActD_trans_table = read.table("~/menghw_HD/our_data/rDNA_data/rDNA_coverage/Hela-n-hic-ActD_w1e5_coverage.bed",header = F,sep = "\t")
# colnames(hela_TR_count) = c("chrom_name","start","end","count","coverage_length","binsize","coverage_rate")
# colnames(hela_IGS_count) = c("chrom_name","start","end","count","coverage_length","binsize","coverage_rate")
# colnames(hela_ActD_TR_count) = c("chrom_name","start","end","count","coverage_length","binsize","coverage_rate")
# colnames(hela_ActD_IGS_count) = c("chrom_name","start","end","count","coverage_length","binsize","coverage_rate")
# colnames(hela_rDNA_trans_table) = c("chrom_name","start","end","count","coverage_length","binsize","coverage_rate")
colnames(hela_rDNA_trans_table) = c("chrom_name","start","end","count","coverage_length","binsize","coverage_rate")
colnames(hela_rDNA_ActD_trans_table) = c("chrom_name","start","end","count","coverage_length","binsize","coverage_rate")

hela_rDNA_hotspot = hela_rDNA_trans_table[hela_rDNA_trans_table$count>=100 & hela_rDNA_trans_table$count <= 2000,]
hela_rDNA_hotspot.ActD = hela_rDNA_ActD_trans_table[hela_rDNA_ActD_trans_table$count>=100 & hela_rDNA_ActD_trans_table$count <= 2000,]

write.table(hela_rDNA_hotspot,file = "~/menghw_HD/data_table/Hela_active_enhancer/Hela_rDNA_hotspot.table",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(hela_rDNA_hotspot.ActD,file = "~/menghw_HD/data_table/Hela_active_enhancer/Hela_rDNA_hotspot_ActD.table",col.names = T,row.names = F,quote = F,sep = "\t")

hela_rDNA_hotspot.merge = read.table(file="~/menghw_HD/data_table/Hela_active_enhancer/Hela_rDNA_hotspot.table.merge",header = F,sep = "\t")
hela_rDNA_hotspot.ActD.merge = read.table(file="~/menghw_HD/data_table/Hela_active_enhancer/Hela_rDNA_hotspot_ActD.table.merge",header = F,sep = "\t")

load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")

hotspot.matrix = matrix(rep(0,46),nrow = 2,ncol = 23)

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom.hotspot.table = hela_rDNA_hotspot.merge[hela_rDNA_hotspot.merge$V1==chrom_name,]
  chrom.hotspot.ActD.table = hela_rDNA_hotspot.ActD.merge[hela_rDNA_hotspot.ActD.merge$V1==chrom_name,]
  hotspot.matrix[1,chrom_index] = sum(chrom.hotspot.table$V3 - chrom.hotspot.table$V2) / chrom.length(chrom_index)
  hotspot.matrix[2,chrom_index] = sum(chrom.hotspot.ActD.table$V3 - chrom.hotspot.ActD.table$V2) / chrom.length(chrom_index)
}


par(mar=c(3,3,2,1),family="Arial",font=1)
barplot(hotspot.matrix,beside = T,col=c("#CC5C55","#42B256"),border = F,space = c(0,0.8),width = 5,xaxt="n",yaxt="n",ylim=c(0,0.3))

axis(side = 2,at=seq(0,0.3,0.1),labels = F,cex.axis=3,lwd=3,tck=-0.05)





