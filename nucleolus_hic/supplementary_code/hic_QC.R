# Hi-C quanlity control code

rm(list=ls())

# data 
{ 
# Hela genome in situ Hi-C
hic_data.CaseName = "Hela-g-hic"
hic_data.RawReadsPair = 4888530936 / 2 / 4
hic_data.MappedReadsPair = hic_data.RawReadsPair * 0.9119
hic_data.MAPQReadsPair = 509042726
hic_data.UndupReadsPair = 212082383 %/% 2 + 347376587
hic_data.RawInter = 347376587
hic_data.ReadsFilter = 319817840
hic_data.FragmentFilter = 247548509 + 55999521
hic_data.Cis = 247548509
hic_data.Trans = 55999521
hic_data.NotSameFragment = 247548509
hic_data.SameFragment = 16269810
}

{ 
  # Hela nucleolus Hi-C
  hic_data.CaseName = "Hela-n-hic"
  hic_data.RawReadsPair = 600e7 / 2 / 4
  hic_data.MappedReadsPair = hic_data.RawReadsPair * 0.9256
  hic_data.MAPQReadsPair = 509042726
  hic_data.UndupReadsPair = 212082383 %/% 2 + 347376587
  hic_data.RawInter = 384757625
  hic_data.ReadsFilter = 342492169
  hic_data.FragmentFilter = 273219336 + 59071557
  hic_data.Cis = 273219336
  hic_data.Trans = 59071557
  hic_data.NotSameFragment = 273219336
  hic_data.SameFragment = 29980981
}


# K562-Dnase-Hic 
{
hic_data.CaseName = "K562-Dnase-hic"
hic_data.RawReadsPair = 1537081408 / 2 / 4
hic_data.MappedReadsPair = hic_data.RawReadsPair * 0.791625 
hic_data.MAPQReadsPair = 236970352 / 2
hic_data.UndupReadsPair = 212035544 / 2
hic_data.RawInter = 78341481
hic_data.ReadsFilter = 78341481
hic_data.FragmentFilter = 78341481
hic_data.Cis = 49605642
hic_data.Trans = 28735839
hic_data.NotSameFragment = 78341481
hic_data.SameFragment = 0
}

# U266-HiC
{
hic_data.CaseName = "U266-hic"
hic_data.RawReadsPair = 834091544 / 2 / 4
hic_data.MappedReadsPair = 208522886 / 2 * 0.9439
hic_data.MAPQReadsPair = hic_data.MappedReadsPair
hic_data.UndupReadsPair = 183274599 / 2
hic_data.RawInter = 85721324
hic_data.ReadsFilter = 77191063
hic_data.FragmentFilter = 43383922 + 23602288
hic_data.Cis = 43383922
hic_data.Trans = 23602288
hic_data.NotSameFragment = 77191063 - 43383922  - 23602288
hic_data.SameFragment = 43383922
}

# U266-FIRRE-C
{
hic_data.CaseName = "U266-F-hic"
hic_data.RawReadsPair = 1471513376 / 2 / 4
hic_data.MappedReadsPair = hic_data.RawReadsPair * 0.8962
hic_data.MAPQReadsPair = hic_data.MappedReadsPair
hic_data.UndupReadsPair = 266032488 / 2
hic_data.RawInter = 101708402
hic_data.ReadsFilter = 95097448
hic_data.FragmentFilter = 53841187 + 29760568
hic_data.Cis = 53841187
hic_data.Trans = 29760568
hic_data.NotSameFragment = 53841187 
hic_data.SameFragment = 65336880 - 53841187
}


# hic_data.RawReadsPair; hic fastq文件的reads pair数目 = fastq文件行数 / 2 / 4
# hic_data.MappedReadsPair; mapped reads数目，可以使用 sam 文件的行数和或者是 hic_data.RawReadsPair * 平均比对率
# hic_data.MAPQReadsPair; 大于MAPQ n 的reads pair的数目，默认是MAPQ20
# hic_data.UndupReadsPair; 去除dup后的reads pair的数目
# hic_data.RawInter; raw cis的interaction的数目
# hic_data.ReadsFilter; 经过reads level过滤后的interaction的数目
# hic_data.FragmentFilter; 经过fragment level过滤后的interaction的数目
# hic_data.Cis; 最终的cis
# hic_data.Trans; 最终的trans
# hic_data.NotSameFragment; 在fragment过滤中，没有在同一个fragment上的interaction
# hic_data.SameFragment; 在fragment过滤中，在同一个fragment上的interaction

# fig1 FastQC read quality distribution
##################################################################################
# fig2 from fastq to raw interaction
##################################################################################
hic_data.vector = c(hic_data.RawReadsPair,
                    hic_data.MappedReadsPair,
                    hic_data.MAPQReadsPair,
                    hic_data.UndupReadsPair,
                    hic_data.RawInter,
                    hic_data.ReadsFilter,
                    hic_data.FragmentFilter
)

## barplot 
barplot.rate = round(hic_data.vector / hic_data.RawReadsPair * 100,2)
barplot.label = paste(barplot.rate,"%",sep = "")

par(mar=c(13,5,5,4),family="Arial",font=2)
barplot(barplot.rate,ylim = c(0,100),space = 1,yaxt="n",col = "#FF8247")

axis(side = 2, at = seq(0,100,length.out =5),labels = paste(seq(0,100,length.out =5),"%",sep = ""),lwd = 3,lwd.ticks = 3,cex.axis=2,tck=-0.05)
text(x = seq(1.5,14,2),y=barplot.rate + 10,srt=0, pos=1, xpd=TRUE,labels=barplot.label,cex=1.5)
text(x = seq(1,13,2),y=par("usr")[3]-2,srt=-60, pos=4, xpd=TRUE,labels=c("Fastq Reads","Mapped Reads","MAPQ Filter","Remove Duplicate","Raw Interaction","Enzyme Filter","Fragment Filter"),cex=2)
title(main = "Hi-C Data Process Rate",cex.main=3)


##################################################################################
# fig3 fragment level filter
##################################################################################
## same fragment v.s. not samefragment
# hic_data.NotSameFragment = 90
# hic_data.SameFragment = 10
barplot.fragment.data = c(hic_data.NotSameFragment,hic_data.SameFragment)
barplot.fragment.rate = round(barplot.fragment.data / sum(barplot.fragment.data),1) * 100
barplot.fragment.label = paste(barplot.fragment.rate,"%",sep = "")
par(mar=c(5,5,5,5))

barplot(barplot.fragment.rate,ylim = c(0,100),yaxt="n",col = c("#FF8247","#7EC0EE"),names.arg = c("Not Same Frag","Same Frag"),cex.names=1.5)
axis(side = 2, at = seq(0,100,length.out =5),labels = paste(seq(0,100,length.out =5),"%",sep = ""),lwd = 3,lwd.ticks = 3,cex.axis=2)
text(x = c(0.5,1.5),y=barplot.fragment.rate + 5,srt=0, pos=4, xpd=TRUE,labels=barplot.fragment.label,cex=2)
title(main = "Hi-C Fragment Filter",cex.main=1.5)


##################################################################################
# fig4 fragment distance distribution
##################################################################################
# U266-g-hic FragDist 
hic_data.FragDist.table = read.table(file="/home/menghw/menghw_HD/data_table/u266_data/U266-g-hic_FragDist_part.txt",header = F)

# U266-f-hic FragDist
hic_data.FragDist.table = read.table(file="/home/menghw/menghw_HD/data_table/u266_data/U266-f-hic_FragDist_part.txt",header = F)

# Hela-n-hic FragDist
hic_data.FragDist.table = read.table(file="/home/menghw/menghw_HD/our_data/Hela-nucleolus-hic-hg19/tmp.data/filtered/Hela-n-hic_FragDist_part.txt",header = F)

# Hela-n-hic FragDist
hic_data.FragDist.table = read.table(file="/home/menghw/menghw_HD/our_data/Hela-genome-hic-hg19/tmp.data/filtered/Hela-g-hic_FragDist.txt",header = F)



hic_data.FragDist.vector = as.vector(hic_data.FragDist.table[,1])
hic_data.FragDist.vector.fix = hic_data.FragDist.vector[hic_data.FragDist.vector>=0]

par(mar=c(5,8,5,5),family="Arial",font=1)
hist(log10(hic_data.FragDist.vector.fix),breaks = seq(0,9,0.1),xaxt="n",yaxt="n",xlab = NULL,ylab = NULL,main = NULL,col = "#7EC0EE",ylim=c(0,35000))

axis(side = 1, at = seq(0,9,2),labels = F,lwd = 3,lwd.ticks = 3,tck=-0.05)
text(x=seq(0,9,2),y=-4000,labels = c("0","0.1K","10K","1M","100M"),cex = 3,xpd=T)

axis(side = 2, at = seq(0,35000,10000),labels = F,lwd = 3,lwd.ticks = 3,cex.axis=2,tck=-0.05)
text(x=-0.5,y=seq(0,35000,10000),labels = seq(0,35000,10000),cex = 3,xpd=T,pos = 2)

text(x=1,y=29500,labels = "Count",cex = 3,xpd=T,srt=0)
box(lwd=3)

title(main = "Fragment Distance Distribution",cex.main=3)

# hic_data.acc.break = c(seq(0,100e3,10e3),seq(100e3,100e6,100e3))
# hic_data.acc.value = hic_data.acc.break
# for(i in c(1:length(hic_data.acc.break))){
#   hic_data.acc.value[i] = length(which(hic_data.FragDist.vector<=hic_data.acc.break[i]))
# }
# hic_data.acc.value.ratio = round(hic_data.acc.value / length(hic_data.FragDist.vector),2)

# par(mar=c(5,5,5,2))
# plot(x=hic_data.acc.break ,hic_data.acc.value.ratio,type="o",xaxt="n",yaxt="n",xlab="Fragment Distance",ylab="Cumulate Rate",cex=1,cex.lab=2,col="#7EC0EE",pch=16,bty="n")
# axis(side = 2, at = seq(0,1,length.out =5),labels = paste(seq(0,100,length.out =5),"%",sep = ""),lwd = 3,lwd.ticks = 3,cex.axis=2)
# axis(side = 1, at = seq(0,100e6,20e6),labels = paste(seq(0,100,20),"Mb",sep = "") ,lwd = 3,lwd.ticks = 3,cex.axis=2)
# box(lwd=3)
# title(main = "Hi-C Fragment Distance Cumulate Plot",cex.main=1.5)


##################################################################################
# fig5 cis trans ratio
##################################################################################
par(mar=c(2,2,5,2))
pie.data = c(hic_data.Cis,hic_data.Trans)
pie.label = c(sprintf("Cis %.1f%%",hic_data.Cis / sum(pie.data) * 100),
              sprintf("Trans %.1f%%",hic_data.Trans / sum(pie.data) * 100))

pie(pie.data,labels = pie.label,col = c("#FF8247","#7EC0EE"),cex=2)
# legend("bottom",c("Intra Reads","Inter Reads"),fill =c("#FF8247","#7EC0EE"),cex=1.2)
title(main = "Cis and Trans Interaction",cex.main=3)

##################################################################################
# fig6 TAD size
##################################################################################

#######################################
## U266-g-hic data
#######################################
rm(list=ls())
binsize = 40000
cutoff = 25

chrom_index = 1
TAD_table.path = sprintf("~/menghw_HD/Project/TAD_calling/u266_TAD_region/U266-g-hic_chr_%d_%d_%1.0f.tad",chrom_index,binsize,cutoff)
TAD_table.all = read.table(TAD_table.path,header = F,sep = "\t")

for(chrom_index in c(2:22)){
  TAD_table.path = sprintf("~/menghw_HD/Project/TAD_calling/u266_TAD_region/U266-g-hic_chr_%d_%d_%1.0f.tad",chrom_index,binsize,cutoff)
  TAD_table.all = rbind(TAD_table.all,read.table(TAD_table.path,header = F,sep = "\t"))
}

TAD_length = (TAD_table.all$V3 - TAD_table.all$V2)
TAD_median = median(TAD_length)
TAD_median = 629e3

# family 是字体，font 1=plain, 2=bold, 3=italic, 4=bold italic, 5=symbol 
par(mar=c(6,10,3,3),family="Arial",font=1)
hist(log10(TAD_length),breaks = seq(4,8,0.2),xlim = c(4,8),ylim=c(0,1200),xlab = "",ylab = "",main="",cex.axis=2,cex.lab=2,lwd=3,col="#804C15",xaxt="n",yaxt="n")

axis(side=1,at=c(4:8),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(4:8),y=-170,labels = c("10K","0.1M","1M","10M","100M"),cex = 3,xpd=T)

axis(side=2,at=seq(0,1200,400),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=3.8,y=seq(0,1200,400),labels = seq(0,1200,400),cex = 3,xpd=T,pos = 2)
text(x=2.8,y=600,labels = "Count",cex = 3,xpd=T,srt=90)

abline(v=log10(TAD_median),lwd=5,pch=10,lty=2,col="black")
text(x=7,y=1100,sprintf("median=%dKb",floor(TAD_median %/% 1000)),pos=1,cex=2.5)
box(lwd=3)



##################################################################################
# fig7 dot plot 两个重复之间的correlation
##################################################################################
rm(list=ls())

###############################################################
# hela-g-hic load
###############################################################
load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-g-hic-rep1_all_matrix_100000.RData")
load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-g-hic-rep2_all_matrix_100000.RData")
load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-g-hic-rep3_all_matrix_100000.RData")

###############################################################
# hela-n-hic load 
###############################################################
load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-n-hic-rep1_all_matrix_100000.RData")
load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-n-hic-rep2_all_matrix_100000.RData")
load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-n-hic-rep3_all_matrix_100000.RData")


###############################################################
# fix seq depth
###############################################################
hela_g_hic_all_matrix.rep1.fix = ceiling(as.matrix(hela_g_hic_all_matrix.rep1 / sum(hela_g_hic_all_matrix.rep1) * 150e6))
hela_g_hic_all_matrix.rep2.fix = ceiling(as.matrix(hela_g_hic_all_matrix.rep2 / sum(hela_g_hic_all_matrix.rep2) * 150e6))
hela_g_hic_all_matrix.rep3.fix = ceiling(as.matrix(hela_g_hic_all_matrix.rep3 / sum(hela_g_hic_all_matrix.rep3) * 150e6))

hela_n_hic_all_matrix.rep1.fix = ceiling(as.matrix(hela_n_hic_all_matrix.rep1 / sum(hela_n_hic_all_matrix.rep1) * 150e6))
hela_n_hic_all_matrix.rep2.fix = ceiling(as.matrix(hela_n_hic_all_matrix.rep2 / sum(hela_n_hic_all_matrix.rep2) * 150e6))
hela_n_hic_all_matrix.rep3.fix = ceiling(as.matrix(hela_n_hic_all_matrix.rep3 / sum(hela_n_hic_all_matrix.rep3) * 150e6))

png(file="~/menghw_HD/R_image/Hela-hic-figure/Hela-g-hic-rep12.png",width = 2000,height = 2000)
par(mar=c(5,5,1,1),font=1)
plot(x=as.vector(hela_g_hic_all_matrix.rep1.fix),y=as.vector(hela_g_hic_all_matrix.rep2.fix),pch=19,col="#5C55DD11",xlim=c(0,300),ylim=c(0,300),xaxt="n",yaxt="n",xlab = "",ylab = "")
axis(side=1,at=seq(0,300,60),labels = F,cex.axis=3,lwd=3,tck=-0.05)
axis(side=2,at=seq(0,300,60),labels = F,cex.axis=3,lwd=3,tck=-0.05)
# text(y=seq(0,200,40),x=-30 ,labels = seq(0,200,40),cex = 2,xpd=T,srt = 0)
abline(a=0,b=1)
box(lwd=3)
dev.off()

cor(as.vector(ceiling(hic_matrix.rep1[1:4000])),as.vector(ceiling(hic_matrix.rep2[1:4000])-3))


png(file="~/menghw_HD/R_image/Hela-hic-figure/Hela-g-hic-rep12.png",width = 2000,height = 2000)
par(mar=c(5,5,1,1),font=1)
plot(x=as.vector(hela_g_hic_all_matrix.rep1.fix),y=as.vector(hela_g_hic_all_matrix.rep2.fix),pch=19,col="#5C55DD11",xlim=c(0,300),ylim=c(0,300),xaxt="n",yaxt="n",xlab = "",ylab = "")
axis(side=1,at=seq(0,300,60),labels = F,cex.axis=3,lwd=3,tck=-0.05)
axis(side=2,at=seq(0,300,60),labels = F,cex.axis=3,lwd=3,tck=-0.05)
# text(y=seq(0,200,40),x=-30 ,labels = seq(0,200,40),cex = 2,xpd=T,srt = 0)
abline(a=0,b=1)
box(lwd=3)
dev.off()

##################################################################################
# fig8 correlation matrix heatmap
##################################################################################
cor_matrix = matrix(rep(0,9),3)






















