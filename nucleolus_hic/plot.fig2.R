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
# hela_n_seq_peak.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_w1e4.table",header = T,sep = "\t")
# hela_n_hic_peak.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix_w1e5.table",header = F,sep = "\t")

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
load("~/menghw_HD/R_code/my_function/MyPlot_v07.RData")
chrom_index = 16
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


##################################################################################
# fig2 boxplot of FPKM
##################################################################################
rm(list=ls())

###########################################
## data loading
###########################################
# gene expression table 
hela_gene_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/raw/Hela_exp_fixed.table",header = T,sep = "\t")

## NAIR table
hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep="\t")
colnames(hela_NAT_table) = c("chrom_name","start","end","name","value")
## NAD table 
hela_NAD_table = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep="\t")
colnames(hela_NAD_table) = c("chrom_name","start","end","name","value")

hela_NAT_complement_table = read.table("~/menghw_HD/data_table/Hela_NAIR_complement.txt",header = F,sep="\t")
hela_NAT_complement_table[,2] = as.numeric(hela_NAT_complement_table[,2])
hela_NAT_complement_table[,3] = as.numeric(hela_NAT_complement_table[,3])
colnames(hela_NAT_complement_table) = c("chrom_name","start","end")

## function
load("~/menghw_HD/R_code/my_function/overlap_gene.RData")

## overlap gene list 
hela_gene_exp_table.fix = hela_gene_exp_table[hela_gene_exp_table$FPKM.control>0.1,]
hela_NAD_gene_exp_table = overlap.gene(region_table = hela_NAD_table,gene_table = hela_gene_exp_table.fix)
hela_NAT_gene_exp_table = overlap.gene(region_table = hela_NAT_table,gene_table = hela_gene_exp_table.fix)
hela_non_NAT_gene_exp_table = overlap.gene(region_table = hela_NAT_complement_table,gene_table = hela_gene_exp_table.fix)

##############################################
## fig2 boxplot FPKM 
#############################################
par(mar=c(5,5,3,1),family="Arial",font=1)
boxplot(log2(hela_NAD_gene_exp_table$FPKM.control),log2(hela_NAT_gene_exp_table$FPKM.control),log2(hela_non_NAT_gene_exp_table$FPKM.control),ylim=c(-5,15),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F,lwd=1.5)
axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
# text(x=c(1:3),y=-8,labels = c("N","E","G"),cex = 3,xpd=T)
# text(x=1.4,y=6.3,labels = "H3K9me3 Count",cex = 3,xpd=T,pos=1)
axis(side=2,at=seq(-5,15,5),labels = F,cex.axis=3,lwd=3,tck=-0.05)
# text(x=0.2,y=seq(-5,15,5),labels =seq(-5,15,5),cex = 3,xpd=T,pos=2)

## statistics test
wilcox.test(log2(hela_NAD_gene_exp_table$FPKM.control),log2(hela_NAT_gene_exp_table$FPKM.control))
wilcox.test(log2(hela_NAT_gene_exp_table$FPKM.control),log2(hela_non_NAT_gene_exp_table$FPKM.control))


##################################################################################
# fig2 boxplot of SINE / LINE&LTR
##################################################################################
chrom_list = c(paste("chr",c(1:22),sep = ""),"chrX")

########################################
## load SINE / LINE LTR coverage table 
########################################
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

###################################
## SINE density
###################################
NAD_SINE_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAD_SINE_coverage.bed",header = F,sep = "\t")
NAT_SINE_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAT_SINE_coverage.bed",header = F,sep = "\t")
# genome_SINE_coverage_table = SINE_coverage_table
other_SINE_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAIR_complement_SINE_coverage.bed",header = F,sep = "\t")

NAD_SINE_coverage_rate = NAD_SINE_coverage_table$V9[NAD_SINE_coverage_table$V9>0]
NAT_SINE_coverage_rate = NAT_SINE_coverage_table$V9[NAT_SINE_coverage_table$V9>0]
other_SINE_coverage_rate = other_SINE_coverage_table$V7
other_SINE_coverage_rate.fix = sample(other_SINE_coverage_rate,size = 2000,replace = F)
# genome_SINE_coverage_rate = genome_SINE_coverage_table$rate[genome_SINE_coverage_table$rate>0]
# genome_SINE_coverage_rate.fix = sample(genome_SINE_coverage_rate,size = 2000,replace = F)


par(mar=c(5,5,3,1),family="Arial",font=1)
boxplot(NAD_SINE_coverage_rate,NAT_SINE_coverage_rate,other_SINE_coverage_rate.fix,ylim=c(0,0.5),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F,lwd=1.5,pch=1)
axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
# text(x=c(1:3),y=-0.06,labels = c("N","E","G"),cex = 3,xpd=T)
# text(x=1.3,y=0.55,labels = "SINE%",cex = 3,xpd=T)
axis(side=2,at=seq(0,0.5,0.25),labels = F,cex.axis=3,lwd=3,tck=-0.05)
# text(x=-0.2,y=seq(0,0.5,0.25),labels = seq(0,50,25),cex = 3,xpd=T)


wilcox.test(NAD_SINE_coverage_rate,NAT_SINE_coverage_rate)
wilcox.test(NAD_SINE_coverage_rate,other_SINE_coverage_rate.fix)


###################################
## LINE/LTR density
###################################
rm(list=ls())
NAD_LINE_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAD_LINE_coverage.bed",header = F,sep = "\t")
NAT_LINE_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAT_LINE_coverage.bed",header = F,sep = "\t")

NAD_LTR_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAD_LTR_coverage.bed",header = F,sep = "\t")
NAT_LTR_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAT_LTR_coverage.bed",header = F,sep = "\t")

other_LINE_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAIR_complement_LINE_coverage.bed",header = F,sep = "\t")
other_LTR_coverage_table = read.table("~/menghw_HD/reference/RepeatMakser/NAIR_complement_LTR_coverage.bed",header = F,sep = "\t")

NAD_LINT_LTR_coverage_rate = NAD_LINE_coverage_table$V9 + NAD_LTR_coverage_table$V9
NAT_LINT_LTR_coverage_rate = NAT_LINE_coverage_table$V9 + NAT_LTR_coverage_table$V9
other_LINE_LTR_coverage_rate = other_LINE_coverage_table$V7 + other_LTR_coverage_table$V7

NAD_LINT_LTR_coverage_rate.fix = NAD_LINT_LTR_coverage_rate
NAT_LINT_LTR_coverage_rate.fix = NAT_LINT_LTR_coverage_rate
other_LINE_LTR_coverage_rate.fix = sample(other_LINE_LTR_coverage_rate,3000,F)


par(mar=c(5,5,3,1),family="Arial",font=1)
boxplot(NAD_LINT_LTR_coverage_rate.fix,NAT_LINT_LTR_coverage_rate.fix,other_LINE_LTR_coverage_rate.fix,ylim=c(0,1),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F,lwd=1.5)

axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
# text(x=c(1:3),y=-0.13,labels = c("N","E","G"),cex = 3,xpd=T)
# text(x=2,y=1.1,labels = "LINE/LTR%",cex = 3,xpd=T,pos = 1)
axis(side=2,at=seq(0,1,0.5),labels = F,cex.axis=3,lwd=3,tck=-0.05)
# text(x=-0.35,y=seq(0,1,0.5),labels = seq(0,100,50),cex = 3,xpd=T)

wilcox.test(NAD_LINT_LTR_coverage_rate.fix,NAT_LINT_LTR_coverage_rate.fix)
wilcox.test(other_LINE_LTR_coverage_rate.fix,NAT_LINT_LTR_coverage_rate.fix)

###################################
## gene density
###################################
NAD_gene_coverage_table = read.table("~/menghw_HD/reference/gene_gtf/NAD_protein_coding_gene_coverage.bed",header = F,sep = "\t")
NAT_gene_coverage_table = read.table("~/menghw_HD/reference/gene_gtf/NAT_protein_coding_gene_coverage.bed",header = F,sep = "\t")
other_gene_coverage_table = read.table("~/menghw_HD/reference/gene_gtf/NAIR_complement_protein_coding_gene_coverage.bed",header = F,sep = "\t")

NAD_gene_density = NAD_gene_coverage_table$V6 / (NAD_gene_coverage_table$V3 - NAD_gene_coverage_table$V2) * 1e6
NAT_gene_density = NAT_gene_coverage_table$V6 / (NAT_gene_coverage_table$V3 - NAT_gene_coverage_table$V2) * 1e6
other_gene_density = other_gene_coverage_table$V4 / (other_gene_coverage_table$V3 - other_gene_coverage_table$V2) * 1e6


wilcox.test(NAD_gene_density[NAD_gene_density>0],NAT_gene_density[NAT_gene_density>0])
wilcox.test(NAD_gene_density[NAD_gene_density>0],genome_gene_density[genome_gene_density>0])

par(mar=c(5,5,3,1),family="Arial",font=1)
boxplot(NAD_gene_density[NAD_gene_density>0],NAT_gene_density[NAT_gene_density>0],other_gene_density[other_gene_density>0],ylim=c(0,60),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F,lwd=1.5)
axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
# text(x=c(1:3),y=-6,labels = c("N","E","G"),cex = 3,xpd=T)
# text(x=1.25,y=56,labels = "No.Gene/Mb",cex = 3,xpd=T)
axis(side=2,at=seq(0,60,20),labels = F,cex.axis=3,lwd=3,tck=-0.05)
# text(x=-0.3,y=seq(0,50,25),labels = seq(0,50,25),cex = 3,xpd=T)
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












