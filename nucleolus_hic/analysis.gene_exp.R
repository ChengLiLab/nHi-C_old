rm(list=ls())
hela_gene_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/raw/Hela_exp_fixed.table",header = T,sep = "\t")
hk_gene_table = read.table("~/menghw_HD/reference/gene_gtf/HK_gene_list.bed",header = F,sep = "\t")
########################################################################
## data loading
########################################################################
## NAT table
hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep="\t")
colnames(hela_NAT_table) = c("chrom_name","start","end","name","value")

hela_NAD_table = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep="\t")
colnames(hela_NAD_table) = c("chrom_name","start","end","name","value")

hela_NAT_complement_table = read.table("~/menghw_HD/data_table/Hela_NAIR_complement.txt",header = F,sep="\t")
hela_NAT_complement_table[,2] = as.numeric(hela_NAT_complement_table[,2])
hela_NAT_complement_table[,3] = as.numeric(hela_NAT_complement_table[,3])
colnames(hela_NAT_complement_table) = c("chrom_name","start","end")
# colnames(hela_NAT_complement_table) = c("chrom_name","start","end","name","value")


## function
# load("~/menghw_HD/R_code/my_function/MyPlot_v08.RData")
load("~/menghw_HD/R_code/my_function/overlap_gene.RData")
########################################################################
## Hela-g-hic analysis
########################################################################
hela_gene_exp_table.fix = hela_gene_exp_table[hela_gene_exp_table$FPKM.control>0,]

hela_NAD_gene_exp_table = overlap.gene(region_table = hela_NAD_table,gene_table = hela_gene_exp_table.fix)
hela_NAT_gene_exp_table = overlap.gene(region_table = hela_NAT_table,gene_table = hela_gene_exp_table.fix)
hela_non_NAT_gene_exp_table = overlap.gene(region_table = hela_NAT_complement_table,gene_table = hela_gene_exp_table.fix)

########################################################################
## fig2 G NAD NAT genome区域表达量的不同
########################################################################
par(mar=c(5,5,3,1),family="Arial",font=1)
boxplot(log2(hela_NAD_gene_exp_table$FPKM.control),log2(hela_NAT_gene_exp_table$FPKM.control),log2(hela_non_NAT_gene_exp_table$FPKM.control),ylim=c(-15,15),col = c("#7080D7","#FF7304","#5CCDC9"),xaxt="n",yaxt="n",frame.plot=F,lwd=1.5)
axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
# text(x=c(1:3),y=-8,labels = c("N","E","G"),cex = 3,xpd=T)
# text(x=1.4,y=6.3,labels = "H3K9me3 Count",cex = 3,xpd=T,pos=1)

axis(side=2,at=seq(-15,15,10),labels = F,cex.axis=3,lwd=3,tck=-0.05)
# text(x=0.2,y=seq(-5,15,5),labels =seq(-5,15,5),cex = 3,xpd=T,pos=2)

wilcox.test(log2(hela_NAD_gene_exp_table$FPKM.control),log2(hela_NAT_gene_exp_table$FPKM.control))
wilcox.test(log2(hela_NAT_gene_exp_table$FPKM.control),log2(hela_non_NAT_gene_exp_table$FPKM.control))




sum(hela_NAT_table$end - hela_NAT_table$start)
sum(hela_NAD_table$end - hela_NAD_table$start)



write.table(hela_NAT_gene_exp_overlap,file = "~/menghw_HD/data_table/Hela-gene_exp_NAT_overlap.bed.v1",col.names = F,row.names = F,quote = F,sep = "\t")
write.table(hela_NAD_gene_exp_overlap,file = "~/menghw_HD/data_table/Hela-gene_exp_NAD_overlap.bed.v1",col.names = F,row.names = F,quote = F,sep = "\t")
write.table(hela_gene_exp_table.fix,file = "~/menghw_HD/data_table/Hela-gene_exp_table.txt",col.names = F,row.names = F,quote = F,sep = "\t")

hela_NAT.hk_gene= overlap(hela_NAT_table,hk_gene_table,level_vetor = c(paste("chr",c(1:22),sep = ""),"chrX"))



sum(hela_NAT_table$end - hela_NAT_table$start) / 2.8e9

NAT_gene_table = NULL
for( gene_name in unique(as.character(hela_NAT_gene_exp_overlap$name.2)) ) {
   NAT_gene_table = rbind(NAT_gene_table,hela_gene_exp_table.fix[hela_gene_exp_table.fix$gene_name==gene_name,])
}

boxplot(NAT_gene_table$fold_change,hela_gene_exp_table.fix$fold_change)

require(vioplot)
vioplot(NAT_gene_table$fold_change,hela_gene_exp_table.fix$fold_change,ylim = c(1:100))


write.table(NAT_gene_table[NAT_gene_table$significant=="yes",],file = "~/menghw_HD/data_table/Hela-gene_exp_NAT_table.bed.v1",col.names = F,row.names = F,quote = F,sep = "\t")


########################################################################
## NAT GO barplot
########################################################################
p_value = c(4.74e-4,4.61e-4,4.84e-4,1.07e-4)
GO_term = c("rRNA processing","Structural constitnent of ribosome","Ribosomal subunit","Mitochodrial inner membrane")
p_value_log = log10(p_value) * (-1)

par(mar=c(8,28,4,2),family="Arial",font=1)
barplot(p_value_log,horiz = T,axes = F,border = F,col = c("#D4AD6A","#5D5393",rep("#44887B",2)),xlim = c(0,4))
axis(side=1,at=seq(0,4),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=seq(0,4),y=-0.7,labels = seq(0,4),cex = 3,xpd=T,pos = NULL)
text(x=2,y=-1.7,labels = "-Log10 P",cex = 3,xpd=T,pos = NULL)

text(x=-0.2,y=c(0.8,1.8,3,4.2),labels = GO_term,cex = 2,xpd=T,srt = 0,pos = 2,offset = 0)

########################################################################
## 表达量
########################################################################

hela_NAT_complement_gene_table = overlap(hela_NAT_complement_table,hela_gene_exp_table.fix,level_vetor = c(paste("chr",c(1:22),sep = ""),"chrX"))

NAT_complement_gene_exp = NULL
for( gene_name in unique(as.character(hela_NAT_complement_gene_table$name.2)) ) {
  NAT_complement_gene_exp = rbind(NAT_complement_gene_exp,hela_gene_exp_table.fix[hela_gene_exp_table.fix$gene_name==gene_name,])
}

boxplot(log10(NAT_gene_table$FPKM.control),log10(NAT_complement_gene_exp$FPKM.control),ylim=c(-1,10))



hist(log10(NAT_complement_gene_exp$FPKM.control))




