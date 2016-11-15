rm(list=ls())

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

SINE_coverage_table.bed.length = SINE_coverage_table[,c(1,2,3,5)]
SINE_coverage_table.bed.density = SINE_coverage_table[,c(1,2,3,7)]

LINE_LTR_coverage_table.bed.length = LINE_LTR_coverage_table[,c(1,2,3,5)]
LINE_LTR_coverage_table.bed.density = LINE_LTR_coverage_table[,c(1,2,3,7)]




chrom.name <- function(chrom_index){
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }else if(chrom_index == 25){
    chrom_name = "rDNA"
  }
}


fix_coordinate <- function(bed_table){
  for(chrom_index in c(1:24)){
    chrom_name = chrom.name(chrom_index)
    chrom_len = chrom.length(chrom_index)
    index = which(bed_table[,1]==chrom_name)[length(which(bed_table[,1]==chrom_name))]  
    bed_table[index,3] = chrom_len - 1
  }
  return(bed_table)
}

chrom.length <- function(chrom_index,GenomeVersion="hg19"){
  # 根据chr_index返回chr_len
  HG19_LEN=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
             135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,
             51304566,155270560,59373566,42999)
  
  HG38_LEN=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,
             135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,
             50818468,156040895,57227415,42999)
  if(GenomeVersion=="hg19"){
    return(HG19_LEN[chrom_index])
  }else if(GenomeVersion=="hg38"){
    return(HG38_LEN[chrom_index])
  }
}

SINE_coverage_table.bed.length = fix_coordinate(SINE_coverage_table.bed.length)
SINE_coverage_table.bed.density = fix_coordinate(SINE_coverage_table.bed.density)

LINE_LTR_coverage_table.bed.length = fix_coordinate(LINE_LTR_coverage_table.bed.length)
LINE_LTR_coverage_table.bed.density = fix_coordinate(LINE_LTR_coverage_table.bed.density)


write.table(SINE_coverage_table.bed.length,file="~/menghw_HD/data_table/WashU_File/hg19_SINE_coverage_length.bedgraph",col.names = F,row.names = F,quote = F,sep = "\t")
write.table(SINE_coverage_table.bed.density,file="~/menghw_HD/data_table/WashU_File/hg19_SINE_coverage_density.bedgraph",col.names = F,row.names = F,quote = F,sep = "\t")

write.table(LINE_LTR_coverage_table.bed.length,file="~/menghw_HD/data_table/WashU_File/hg19_LINE_LTR_coverage_length.bedgraph",col.names = F,row.names = F,quote = F,sep = "\t")
write.table(LINE_LTR_coverage_table.bed.density,file="~/menghw_HD/data_table/WashU_File/hg19_LINE_LTR_coverage_density.bedgraph",col.names = F,row.names = F,quote = F,sep = "\t")

####################################################################################
# Fig 2A plot
####################################################################################
load(file="~/menghw_HD/R_code/my_function/MyPlot_v05.RData")
chrom_index = 16
chrom_start = 0
chrom_end = chrom.length(chrom_index)

hela_NAT = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
hela_NAD = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep = "\t")

SINE_coverage_table.bed.length.fix = SINE_coverage_table.bed.length
SINE_coverage_table.bed.length.fix$cover_len = SINE_coverage_table.bed.length$cover_len - 10000
SINE_coverage_table.bed.length.fix$cover_len[SINE_coverage_table.bed.length.fix$cover_len<0] = 0

LINE_LTR_coverage_table.bed.length.fix = LINE_LTR_coverage_table.bed.length
LINE_LTR_coverage_table.bed.length.fix$cover_len = LINE_LTR_coverage_table.bed.length$cover_len - 20000
LINE_LTR_coverage_table.bed.length.fix$cover_len[LINE_LTR_coverage_table.bed.length.fix$cover_len<0] = 0


fig.facet <- layout(matrix(c(1:4),nrow = 4,byrow = FALSE),width = rep(10,4),heights = rep(1,4))
layout.show(fig.facet)

par(mar=c(0,5,0,1))
plot.barplot(SINE_coverage_table.bed.length.fix,chrom_index,chrom_start,chrom_end,ylim = c(0,6e4),track_color = "#B2200E")
axis(side=2,at=c(0,5e4),labels = F,lwd=3)
text(x=-6e6,y=c(3000,4.7e4) ,labels =c("0","5E4") ,cex = 2,xpd=T,srt = 0)
text(x=-3e6,y=4.5e4 ,labels ="Coverage" ,cex = 2,xpd=T,srt = 0,pos=4)

# box(lwd=2)
par(new=T)
plot.peak(df.peak = hela_NAT,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF730455",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)


par(mar=c(0,5,0,1))
plot.barplot(LINE_LTR_coverage_table.bed.length.fix,chrom_index,chrom_start,chrom_end,ylim = c(0,6e4),track_color = "#805E15")
axis(side=2,at=c(0,5e4),labels = F,lwd=3)
text(x=-6e6,y=c(3000,4.7e4) ,labels =c("0","5E4") ,cex = 2,xpd=T,srt = 0)
text(x=-3e6,y=4.5e4 ,labels ="Coverage" ,cex = 2,xpd=T,srt = 0,pos=4)

par(new=T)
plot.peak(df.peak = hela_NAT,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF730455",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
# box(lwd=2)


hg19_RNA_gene_cover_table = read.table("~/menghw_HD/reference/gene_gtf/hg19_pro_coding_coverage_100000.bed",header = F,sep = "\t")
hg19_gene_table = hg19_RNA_gene_cover_table[,c(1,2,3,4)]
hg19_gene_table.fix = hg19_gene_table
# hg19_gene_table.fix$V5 = hg19_gene_table$V5 - median(hg19_gene_table$V5)
# hg19_gene_table.fix$V5[hg19_gene_table.fix$V5 <0 ] = 0
# quantile(hg19_gene_table.fix$V5)
hg19_gene_table.fix$V4 = hg19_gene_table$V4
par(mar=c(0,5,0,1))
plot.barplot(hg19_gene_table.fix,chrom_index,chrom_start,chrom_end,ylim = c(0,15),track_color = "#012E34")
# box(lwd=2)
axis(side=2,at=c(0,13),labels = F,lwd=3)
text(x=-6e6,y=c(2,11) ,labels =c("0","15") ,cex = 2,xpd=T,srt = 0)
text(x=-3e6,y=11 ,labels ="Count" ,cex = 2,xpd=T,srt = 0,pos=4)
par(new=T)
plot.peak(df.peak = hela_NAT,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF730455",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

# RNA gene table 趋势不明显！！！！
hg19_RNA_gene_table = read.table("~/menghw_HD/reference/gene_gtf/hg19_RNA_gene_coverage_100000.bed",header = F,sep = "\t")
hg19_gene_table = hg19_RNA_gene_table[,c(1,2,3,5)]
hg19_gene_table.fix = hg19_gene_table
# hg19_gene_table.fix$V5 = hg19_gene_table$V5 - median(hg19_gene_table$V5)
# hg19_gene_table.fix$V5[hg19_gene_table.fix$V5 <0 ] = 0
# quantile(hg19_gene_table.fix$V5)
hg19_gene_table.fix$V4 = hg19_gene_table$V4
par(mar=c(0,5,0,1))
plot.barplot(hg19_gene_table.fix,chrom_index,chrom_start,chrom_end,ylim = c(0,100e3),track_color = "#0E464E")
# box(lwd=2)
axis(side=2,at=c(0,13),labels = F,lwd=3)
text(x=-6e6,y=c(2,11) ,labels =c("0","15") ,cex = 2,xpd=T,srt = 0)
text(x=-3e6,y=11 ,labels ="Count" ,cex = 2,xpd=T,srt = 0,pos=4)
par(new=T)
plot.peak(df.peak = hela_NAT,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF730488",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

# Histone Data coverage 
#### active
hela_H3K4me3_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_H3K4me3_ENCFF000BCO_coverage_100000.bed",header = F,sep = "\t")
# hela_H3K36me3_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_H3K36me3_ENCFF000BCA_coverage_100000.bed",header = F,sep = "\t")
# hela_POL2_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_POL2_ENCFF000PEL_coverage_100000.bed",header = F,sep = "\t")
#### repressive
hela_H3K27me3_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_H3K27me3_ENCFF000BBS_coverage_100000.bed",header = F,sep = "\t")
hela_H3K9me3_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_H3K9me3_ENCFF000BBG_coverage_100000.bed",header = F,sep = "\t")
# hela_H4K20me3_table = read.table("~/menghw_HD/data_table/public_table/Hela-ChIP_H4K20me1_ENCFF000BDC_coverage_100000.bed",header = F,sep = "\t")

# active signal H3K4me3
hela_H3K4me3_table.fix = hela_H3K4me3_table[,c(1,2,3,4)]
# hela_H3K4me3_table.fix$V5 = hela_H3K4me3_table.fix$V5 - median(hela_H3K4me3_table$V5)
# hela_H3K4me3_table.fix$V5[hela_H3K4me3_table.fix$V5<0] = 0
hela_H3K4me3_table.fix$V4 = hela_H3K4me3_table.fix$V4 - median(hela_H3K4me3_table$V4)
hela_H3K4me3_table.fix$V4[hela_H3K4me3_table.fix$V4<0] = 0
hela_H3K4me3_table.fix$V4[hela_H3K4me3_table.fix$V4>1800] = 1800
par(mar=c(0,5,0,1))
plot.barplot(hela_H3K4me3_table.fix,chrom_index,chrom_start,chrom_end,ylim = c(0,2000),track_color = "#CC7526")
# box(lwd=2)
axis(side=2,at=c(0,1800),labels = F,lwd=3)
text(x=-6e6,y=c(200,1600) ,labels =c("0","2E3") ,cex = 2,xpd=T,srt = 0)
text(x=-3e6,y=1600 ,labels ="Count" ,cex = 2,xpd=T,srt = 0,pos=4)

par(new=T)
plot.peak(df.peak = hela_NAT,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF730455",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

# repressive signal H3K9me3
hela_H3K9me3_table.fix = hela_H3K9me3_table[,c(1,2,3,4)]
# hela_H3K27me3_table.fix$V5 = hela_H3K27me3_table.fix$V5 - median(hela_H3K27me3_table.fix$V5)
# hela_H3K27me3_table.fix$V5[hela_H3K27me3_table.fix$V5<0] = 0
hela_H3K9me3_table.fix$V4 = hela_H3K9me3_table.fix$V4 - median(hela_H3K9me3_table.fix$V4)
hela_H3K9me3_table.fix$V4[hela_H3K9me3_table.fix$V4<0] = 0
hela_H3K9me3_table.fix$V4[hela_H3K9me3_table.fix$V4>1600] = 1600

par(mar=c(0,5,0,1))
plot.barplot(hela_H3K9me3_table.fix,chrom_index,chrom_start,chrom_end,ylim = c(0,2e3),track_color = "#5820FF")
# box(lwd=2)
axis(side=2,at=c(0,1800),labels = F,lwd=3)
text(x=-6e6,y=c(200,1600) ,labels =c("0","2E3") ,cex = 2,xpd=T,srt = 0)
text(x=-3e6,y=1600 ,labels ="Count" ,cex = 2,xpd=T,srt = 0,pos=4)

par(new=T)
plot.peak(df.peak = hela_NAT,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF730455",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

load(file="~/menghw_HD/R_code/my_function/MyPlot_v05.RData")
par(mar=c(0,5,0,1))
plot.chromosome(chrom_index = 16,chrom_start = 0,chrom_end = chrom_end,border.lwd = 3)

####################################################################################
# Fig 2 B,C,D,E,F plot
####################################################################################
chrom_list = c(paste("chr",c(1:22),sep = ""),"chrX")

###################################
## other region
###################################


###################################
## SINE density
###################################
NAD_SINE_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAD_SINE_coverage.bed",header = F,sep = "\t")
NAT_SINE_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAT_SINE_coverage.bed",header = F,sep = "\t")
genome_SINE_coverage_table = SINE_coverage_table

NAD_SINE_coverage_rate = NAD_SINE_coverage_table$V9[NAD_SINE_coverage_table$V9>0]
NAT_SINE_coverage_rate = NAT_SINE_coverage_table$V9[NAT_SINE_coverage_table$V9>0]
genome_SINE_coverage_rate = genome_SINE_coverage_table$rate[genome_SINE_coverage_table$rate>0]
genome_SINE_coverage_rate.fix = sample(genome_SINE_coverage_rate,size = 2000,replace = F)

par(mar=c(5,5,3,1),family="Arial",font=1)
boxplot(NAD_SINE_coverage_rate,NAT_SINE_coverage_rate,genome_SINE_coverage_rate.fix,ylim=c(0,0.5),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F,lwd=1.5,pch=1)

axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(1:3),y=-0.06,labels = c("N","E","G"),cex = 3,xpd=T)

# text(x=1.3,y=0.55,labels = "SINE%",cex = 3,xpd=T)
axis(side=2,at=seq(0,0.5,0.25),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=-0.2,y=seq(0,0.5,0.25),labels = seq(0,50,25),cex = 3,xpd=T)


wilcox.test(NAD_SINE_coverage_rate,NAT_SINE_coverage_rate)
wilcox.test(NAD_SINE_coverage_rate,genome_SINE_coverage_rate.fix)


###################################
## LINE/LTR density
###################################
NAD_LINE_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAD_LINE_coverage.bed",header = F,sep = "\t")
NAT_LINE_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAT_LINE_coverage.bed",header = F,sep = "\t")

NAD_LTR_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAD_LTR_coverage.bed",header = F,sep = "\t")
NAT_LTR_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAT_LTR_coverage.bed",header = F,sep = "\t")

genome_LINE_LTR_coverage_table = LINE_LTR_coverage_table

NAD_LINT_LTR_coverage_rate = NAD_LINE_coverage_table$V9 + NAD_LTR_coverage_table$V9
NAT_LINT_LTR_coverage_rate = NAT_LINE_coverage_table$V9 + NAT_LTR_coverage_table$V9
genome_LINE_LTR_coverage_rate = genome_LINE_LTR_coverage_table$rate


NAD_LINT_LTR_coverage_rate.fix = NAD_LINT_LTR_coverage_rate[NAD_LINT_LTR_coverage_rate>0]
NAT_LINT_LTR_coverage_rate.fix = NAT_LINT_LTR_coverage_rate[NAT_LINT_LTR_coverage_rate>0]
genome_LINE_LTR_coverage_rate.fix = sample(genome_LINE_LTR_coverage_rate[genome_LINE_LTR_coverage_rate>0],2000,F)


par(mar=c(5,5,3,1),family="Arial",font=1)
boxplot(NAD_LINT_LTR_coverage_rate.fix,NAT_LINT_LTR_coverage_rate.fix,genome_LINE_LTR_coverage_rate.fix,ylim=c(0,1),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F,lwd=1.5)

axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(1:3),y=-0.13,labels = c("N","E","G"),cex = 3,xpd=T)

# text(x=2,y=1.1,labels = "LINE/LTR%",cex = 3,xpd=T,pos = 1)

axis(side=2,at=seq(0,1,0.5),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=-0.35,y=seq(0,1,0.5),labels = seq(0,100,50),cex = 3,xpd=T)

wilcox.test(NAD_LINT_LTR_coverage_rate.fix,NAT_LINT_LTR_coverage_rate.fix)
wilcox.test(genome_LINE_LTR_coverage_rate.fix,NAT_LINT_LTR_coverage_rate.fix)

###################################
## gene density
###################################
NAD_gene_coverage_table = read.table("~/menghw_HD/reference/gene_gtf/NAD_protein_coding_gene_coverage.bed",header = F,sep = "\t")
NAT_gene_coverage_table = read.table("~/menghw_HD/reference/gene_gtf/NAT_protein_coding_gene_coverage.bed",header = F,sep = "\t")
hg19_gene_coverage_table = read.table("~/menghw_HD/reference/gene_gtf/hg19_pro_coding_coverage_100000.bed",header = F,sep = "\t")

NAD_gene_density = NAD_gene_coverage_table$V6 / (NAD_gene_coverage_table$V3 - NAD_gene_coverage_table$V2) * 1e6
NAT_gene_density = NAT_gene_coverage_table$V6 / (NAT_gene_coverage_table$V3 - NAT_gene_coverage_table$V2) * 1e6
genome_gene_density = hg19_gene_coverage_table$V4 / (hg19_gene_coverage_table$V3 - hg19_gene_coverage_table$V2) * 1e6

genome_gene_density.fix = sample(x = genome_gene_density[genome_gene_density>0.1],size = 1000,replace = T)


wilcox.test(NAD_gene_density[NAD_gene_density>0],NAT_gene_density[NAT_gene_density>0])
wilcox.test(NAD_gene_density[NAD_gene_density>0],genome_gene_density[genome_gene_density>0])

par(mar=c(5,5,3,1),family="Arial",font=1)
boxplot(NAD_gene_density[NAD_gene_density>0],NAT_gene_density[NAT_gene_density>0],genome_gene_density.fix,ylim=c(0,50),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F,lwd=1.5)
axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(1:3),y=-6,labels = c("N","E","G"),cex = 3,xpd=T)
# text(x=1.25,y=56,labels = "No.Gene/Mb",cex = 3,xpd=T)

axis(side=2,at=seq(0,50,25),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=-0.3,y=seq(0,50,25),labels = seq(0,50,25),cex = 3,xpd=T)
# box(lwd=3)


###################################
## H3K4me3 signal
###################################
NAD_histone_table = read.table("~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/NAD_H3K4me3_coverage.bed",header = F,sep = "\t")
NAT_histone_table = read.table("~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/NAT_H3K4me3_coverage.bed",header = F,sep = "\t")
genome_histone_table = hela_H3K4me3_table

NAD_histone_density = NAD_histone_table$V6 / (NAD_histone_table$V3 - NAD_histone_table$V2) * 10e6
NAT_histone_density = NAT_histone_table$V6 / (NAT_histone_table$V3 - NAT_histone_table$V2) * 10e6 
genome_histone_density = genome_histone_table$V4 / (genome_histone_table$V3 - genome_histone_table$V2) * 10e6 

NAD_histone_density.fix = log10(NAD_histone_density[NAD_histone_density>0])
NAT_histone_density.fix = log10(NAT_histone_density[NAT_histone_density>0])
genome_histone_density.fix = log10(sample(genome_histone_density[genome_histone_density>0],5000))


par(mar=c(5,5,3,1),family="Arial",font=1)
boxplot(NAD_histone_density.fix,NAT_histone_density.fix,genome_histone_density.fix,ylim=c(4,6),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F,lwd=1.5)
axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)

text(x=c(1:3),y=3.75,labels = c("N","E","G"),cex = 3,xpd=T)
# text(x=1.4,y=6.3,labels = "H3K4me3 Count",cex = 3,xpd=T,pos=1)

axis(side=2,at=seq(4,6,1),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=0,y=seq(4,6,1),labels =c("4","5","6"),cex = 3,xpd=T)

wilcox.test(NAD_histone_density.fix,NAT_histone_density.fix)
wilcox.test(genome_histone_density.fix,NAT_histone_density.fix)

###################################
## H3K9me3 signal
###################################
NAD_histone_table = read.table("~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/NAD_H3K9me3_coverage.bed",header = F,sep = "\t")
NAT_histone_table = read.table("~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/NAT_H3K9me3_coverage.bed",header = F,sep = "\t")
genome_histone_table = hela_H3K9me3_table

NAD_histone_density = NAD_histone_table$V6 / (NAD_histone_table$V3 - NAD_histone_table$V2) * 10e6
NAT_histone_density = NAT_histone_table$V6 / (NAT_histone_table$V3 - NAT_histone_table$V2) * 10e6 
genome_histone_density = genome_histone_table$V4 / (genome_histone_table$V3 - genome_histone_table$V2) * 10e6 /1.1

NAD_histone_density.fix = log10(NAD_histone_density[NAD_histone_density>0])
NAT_histone_density.fix = log10(NAT_histone_density[NAT_histone_density>0])
genome_histone_density.fix = log10(sample(genome_histone_density[genome_histone_density>0],3000))


par(mar=c(5,5,3,1),family="Arial",font=1)
boxplot(NAD_histone_density.fix,NAT_histone_density.fix,genome_histone_density.fix,ylim=c(4,6),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F,lwd=1.5)
axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)

text(x=c(1:3),y=3.75,labels = c("N","E","G"),cex = 3,xpd=T)
# text(x=1.4,y=6.3,labels = "H3K9me3 Count",cex = 3,xpd=T,pos=1)

axis(side=2,at=seq(4,6,1),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=0,y=seq(4,6,1),labels =c(4:6),cex = 3,xpd=T)

wilcox.test(NAD_histone_density.fix,NAT_histone_density.fix)
wilcox.test(genome_histone_density.fix,NAT_histone_density.fix)

