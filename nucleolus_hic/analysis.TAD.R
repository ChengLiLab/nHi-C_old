#######################################
## U266-f-hic data
#######################################
rm(list=ls())
binsize = 40000
cutoff = 25

chrom_index = 1
TAD_table.path = sprintf("~/menghw_HD/Project/TAD_calling/u266_TAD_region/U266-f-hic_chr_%d_%d_%1.0f.tad",chrom_index,binsize,cutoff)
TAD_table.all = read.table(TAD_table.path,header = F,sep = "\t")

for(chrom_index in c(2:22)){
  TAD_table.path = sprintf("~/menghw_HD/Project/TAD_calling/u266_TAD_region/U266-f-hic_chr_%d_%d_%1.0f.tad",chrom_index,binsize,cutoff)
  TAD_table.all = rbind(TAD_table.all,read.table(TAD_table.path,header = F,sep = "\t"))
}

TAD_length = (TAD_table.all$V3 - TAD_table.all$V2)
TAD_median = median(TAD_length)
TAD_median = 320*1e3

# family 是字体，font 1=plain, 2=bold, 3=italic, 4=bold italic, 5=symbol 
par(mar=c(6,10,3,3),family="Arial",font=2)
hist(log10(TAD_length),breaks = seq(4,8,0.2),xlim = c(4,8),ylim=c(0,1200),xlab = "",ylab = "",main="",cex.axis=2,cex.lab=2,lwd=3,col="#FF7304",xaxt="n",yaxt="n")

axis(side=1,at=c(4:8),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(4:8),y=-170,labels = c("10K","0.1M","1M","10M","100M"),cex = 3,xpd=T)

axis(side=2,at=seq(0,1200,400),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=3.8,y=seq(0,1200,400),labels = seq(0,1200,400),cex = 3,xpd=T,pos = 2)
text(x=2.8,y=600,labels = "Count",cex = 3,xpd=T,srt=90)

abline(v=log10(TAD_median),lwd=5,pch=10,lty=2,col="black")
text(x=7,y=1100,sprintf("median=%dKb",ceiling(TAD_median %/% 1000)),pos=1,cex=2.5)
box(lwd=3)


tad_table = read.table(file = "~/menghw_HD/Project/TAD_calling/test.tad",header = F,sep = "\t")
quantile(tad_table$V3 - tad_table$V2)

hela_n_hic_tad_table = read.table(file = "~/menghw_HD/Project/TAD_calling/test.tad.2",header = F,sep = "\t")
quantile(hela_n_hic_tad_table$V3 - hela_n_hic_tad_table$V2)

hela_n_hic_ActD_tad_table = read.table(file = "~/menghw_HD/Project/TAD_calling/test.tad.3",header = F,sep = "\t")
quantile(hela_n_hic_ActD_tad_table$V3 - hela_n_hic_ActD_tad_table$V2)


#######################################
## Hela-n-hic-ActD data
#######################################
rm(list=ls())
binsize = 20000
cutoff = 65

chrom_index = 1
TAD_table.path = sprintf("~/menghw_HD/Project/TAD_calling/TAD_region/Hela-n-hic-ActD_chr_%d_%d_%1.0f.tad",chrom_index,binsize,cutoff)
TAD_table.all = read.table(TAD_table.path,header = F,sep = "\t")

for(chrom_index in c(2:22)){
  TAD_table.path = sprintf("~/menghw_HD/Project/TAD_calling/TAD_region/Hela-n-hic-ActD_chr_%d_%d_%1.0f.tad",chrom_index,binsize,cutoff)
  TAD_table.all = rbind(TAD_table.all,read.table(TAD_table.path,header = F,sep = "\t"))
}

chrom_index = 1
TAD_table.path = sprintf("~/menghw_HD/Project/TAD_calling/TAD_region/Hela-g-hic_chr_%d_%d_%1.0f.tad",chrom_index,binsize,cutoff)
TAD_table.all = read.table(TAD_table.path,header = F,sep = "\t")

for(chrom_index in c(2:22)){
  TAD_table.path = sprintf("~/menghw_HD/Project/TAD_calling/TAD_region/Hela-g-hic_chr_%d_%d_%1.0f.tad",chrom_index,binsize,cutoff)
  TAD_table.all = rbind(TAD_table.all,read.table(TAD_table.path,header = F,sep = "\t"))
}

TAD_length = (TAD_table.all$V3 - TAD_table.all$V2)
TAD_median = median(TAD_length)
TAD_median = 320*1e3

# family 是字体，font 1=plain, 2=bold, 3=italic, 4=bold italic, 5=symbol 
par(mar=c(6,10,3,3),family="Arial",font=2)
hist(log10(TAD_length),breaks = seq(4,8,0.2),xlim = c(4,8),ylim=c(0,1200),xlab = "",ylab = "",main="",cex.axis=2,cex.lab=2,lwd=3,col="#FF7304",xaxt="n",yaxt="n")

axis(side=1,at=c(4:8),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(4:8),y=-170,labels = c("10K","0.1M","1M","10M","100M"),cex = 3,xpd=T)

axis(side=2,at=seq(0,1200,400),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=3.8,y=seq(0,1200,400),labels = seq(0,1200,400),cex = 3,xpd=T,pos = 2)
text(x=2.8,y=600,labels = "Count",cex = 3,xpd=T,srt=90)

abline(v=log10(TAD_median),lwd=5,pch=10,lty=2,col="black")
text(x=7,y=1100,sprintf("median=%dKb",ceiling(TAD_median %/% 1000)),pos=1,cex=2.5)
box(lwd=3)


tad_table = read.table(file = "~/menghw_HD/Project/TAD_calling/test.tad",header = F,sep = "\t")
quantile(tad_table$V3 - tad_table$V2)

hela_n_hic_tad_table = read.table(file = "~/menghw_HD/Project/TAD_calling/test.tad.2",header = F,sep = "\t")
quantile(hela_n_hic_tad_table$V3 - hela_n_hic_tad_table$V2)

hela_n_hic_ActD_tad_table = read.table(file = "~/menghw_HD/Project/TAD_calling/test.tad.3",header = F,sep = "\t")
quantile(hela_n_hic_ActD_tad_table$V3 - hela_n_hic_ActD_tad_table$V2)


##############################################################################
## NAIR and TAD overlap
##############################################################################
rm(list=ls())

## TAD table
g_hic_TAD = NULL
for(i in c(1:23)){
  TAD_file.name = sprintf("~/menghw_HD/Project/TAD_calling/TAD_region/Hela-g-hic_chr_%d_20000_65.tad",i)
  TAD_table = read.table(TAD_file.name,header = F,sep = "\t")
  g_hic_TAD = rbind(g_hic_TAD,TAD_table)
}

g_hic_TAD = cbind(g_hic_TAD,rep("TAD",nrow(g_hic_TAD)))
g_hic_TAD = cbind(g_hic_TAD,rep(1,nrow(g_hic_TAD)))
colnames(g_hic_TAD) = c("chrom_name","start","end","name","value")

# write.table(g_hic_TAD,"~/menghw_HD/data_table/TAD_region/Hela-n-hic_20000_TAD_norm.bed",col.names = F,row.names = F,quote = F,sep = "\t")

## NAIR table
hela_NAIR_table = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix_w1e5.table",header = F,sep="\t")
colnames(hela_NAIR_table) = c("chrom_name","start","end","name","value")

## load function
load(file = "~/menghw_HD/R_code/my_function/MyPlot_v03.RData")

NAIR_TAD_overlap = overlap(hela_NAIR_table,g_hic_TAD,level_vetor = c(paste("chr",c(1:22),sep = ""),"chrX"))
# save.image(file = "~/menghw_HD/R_image/Hela-TAD-NAIR.RData")

sum(NAIR_TAD_overlap$end - NAIR_TAD_overlap$start) / sum(hela_NAIR_table$end - hela_NAIR_table$start)

sum(as.numeric(NAIR_TAD_overlap$end.2) - as.numeric(NAIR_TAD_overlap$start.2))
sum(as.numeric(g_hic_TAD$end) - as.numeric(g_hic_TAD$start))

write.table(g_hic_TAD,file = "~/menghw_HD/data_table/Hela-g-hic_TAD_raw.bed",col.names = T,row.names = F,quote = F,sep = "\t")

##############################################################################
## plot NAIR and TAD 
##############################################################################
load("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_14_20000.RData")

chrom_index = 14
chrom_start = 25e6
chrom_end = 30e6
binsize = 20e3

hela_g_hic_inter = convert_matrix2inter.cis(g_hic_matrix.norm,chrom_index,chrom_start,chrom_end,binsize)

hic_inter_table = hela_g_hic_inter
log_scale = T
WashU_pairwise.fmt = "%s:%d,%d\t%s:%d,%d\t%.4f"
file_connect = file("~/menghw_HD/data_table/WashU_File/Hela-g-hic_chr14_200000_25e6_30e6.norm.log",open = "w")

for(i in c(1:nrow(hic_inter_table))){
  if(i %% 10000==0){
    print(i)  
  }
  inter_table.row = hic_inter_table[i,]
  inter_table.str = NULL
  if(log_scale){
    inter_table.str = sprintf(WashU_pairwise.fmt,inter_table.row[1,1],inter_table.row[1,2],inter_table.row[1,3],inter_table.row[1,4],inter_table.row[1,5],inter_table.row[1,6],inter_table.row[1,8])
  }else{
    inter_table.str = sprintf(WashU_pairwise.fmt,inter_table.row[1,1],inter_table.row[1,2],inter_table.row[1,3],inter_table.row[1,4],inter_table.row[1,5],inter_table.row[1,6],inter_table.row[1,7])
  }
  writeLines(inter_table.str,con = file_connect)
}
close(file_connect)


hela_n_hic_inter = convert_matrix2inter.cis(n_hic_matrix,chrom_index,chrom_start,chrom_end,binsize)

hic_inter_table = hela_n_hic_inter
log_scale = F
WashU_pairwise.fmt = "%s:%d,%d\t%s:%d,%d\t%.4f"
file_connect = file("~/menghw_HD/data_table/WashU_File/Hela-n-hic_chr14_20000_25e6_30e6.raw",open = "w")

for(i in c(1:nrow(hic_inter_table))){
  if(i %% 10000==0){
    print(i)  
  }
  inter_table.row = hic_inter_table[i,]
  inter_table.str = NULL
  if(log_scale){
    inter_table.str = sprintf(WashU_pairwise.fmt,inter_table.row[1,1],inter_table.row[1,2],inter_table.row[1,3],inter_table.row[1,4],inter_table.row[1,5],inter_table.row[1,6],inter_table.row[1,8])
  }else{
    inter_table.str = sprintf(WashU_pairwise.fmt,inter_table.row[1,1],inter_table.row[1,2],inter_table.row[1,3],inter_table.row[1,4],inter_table.row[1,5],inter_table.row[1,6],inter_table.row[1,7])
  }
  writeLines(inter_table.str,con = file_connect)
}
close(file_connect)

##############################################################################
## TAD length distribution
##############################################################################
colnames(hela_g_hic_TAD) = c("chrom_name","start","end","insulation","value")
hela_g_hic_TAD$start = as.numeric(hela_g_hic_TAD$start)
hela_g_hic_TAD$end = as.numeric(hela_g_hic_TAD$end)
hela_g_hic_TAD.fix = hela_g_hic_TAD[(hela_g_hic_TAD$end-hela_g_hic_TAD$start)>300e3,]


sum(hela_g_hic_TAD.fix$end - hela_g_hic_TAD.fix$start)
quantile((hela_g_hic_TAD.fix$end - hela_g_hic_TAD.fix$start))
hist(log10(hela_g_hic_TAD.fix$end - hela_g_hic_TAD.fix$start))

dev.off()

TAD_length = (hela_g_hic_TAD.fix$end - hela_g_hic_TAD.fix$start)
TAD_median = median((hela_g_hic_TAD.fix$end - hela_g_hic_TAD.fix$start))

par(mar=c(6,10,3,3),family="Arial",font=1)
hist(log10(TAD_length),breaks = seq(5,8,0.2),xlim = c(5,7),ylim=c(0,1200),xlab = "",ylab = "",main="",cex.axis=2,cex.lab=2,lwd=3,col="#00685A",xaxt="n",yaxt="n")

axis(side=1,at=c(5:7),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
text(x=c(5:7),y=-200,labels = c("100Kb","1M","10Mb"),cex = 3,xpd=T)

axis(side=2,at=seq(0,1200,400),labels = F,cex.axis=3,lwd=3,tck=-0.05)
text(x=4.8,y=seq(0,1200,400),labels = seq(0,1200,400),cex = 3,xpd=T,pos = 2)
text(x=3.8,y=600,labels = "Count",cex = 3,xpd=T,srt=90)

abline(v=log10(TAD_median),lwd=5,pch=10,lty=2,col="black")
text(x=6.5,y=1100,sprintf("median=%dKb",ceiling(TAD_median %/% 1000)),pos=1,cex=2.5)

box(lwd=3)




##############################################################################
# Fig3 A
##############################################################################
load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")

chrom_index = 16
binsize = 20e3
image.path = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_%d_%d.RData",chrom_index,binsize)
load(file=image.path)

hela_NAT = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
hela_NAD = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep = "\t")

hela_CTCF = read.table("~/menghw_HD/public_data/Hela-ChIP-Seq-data/raw_data/Hela-ChIP_CTCF_ENCFF000BAJ_coverage_10000.bed",header = F,sep = '\t')
hela_CTCF.fix = hela_CTCF[,c(1,2,3,5)]
hela_CTCF.fix$V5 = hela_CTCF$V5 - median(hela_CTCF$V5)
hela_CTCF.fix = hela_CTCF.fix[hela_CTCF.fix$V5>=0,]
# hela_g_hic_TAD = read.table("~/menghw_HD/data_table/TAD_region/Hela-g-hic_TAD_region.bed",header = F,sep = "\t")
# hela_g_hic_TAD=cbind(hela_g_hic_TAD,rep(1,nrow(hela_g_hic_TAD)))

# hela_g_hic_TAD.overlap.raw = overlap(hela_NAT,hela_g_hic_TAD,level_vetor = c(paste("chr",c(1:22),sep = ""),"chrX"))
# hela_g_hic_TAD.overlap.raw.fix = hela_g_hic_TAD.overlap.raw[,c(1,8,9,4,10)]


hela_NAT.fix = hela_NAT
# hela_NAT.fix$V2 = (hela_NAT$V2 %/% 100e3) * 100e3 
# hela_NAT.fix$V3 = (hela_NAT$V3 %/% 100e3) * 100e3 + 50e3


chrom_start = 22e6
chrom_end = 28e6
g_hic_matrix.part = matrix.part(g_hic_matrix.norm,chrom_index,chrom_start,chrom_end,binsize)
n_hic_matrix.part = matrix.part(n_hic_matrix,chrom_index,chrom_start,chrom_end,binsize)

fig.facet <- layout(matrix(c(1:4),nrow = 4,byrow = FALSE),width = rep(10,4),heights = c(7,1,7,2))
layout.show(fig.facet)

par(mar=c(0.5,2,0.5,2))
plot.TAD(g_hic_matrix.part,ylim=c(0,100),maxBound = 0.97)

par(mar=c(0.5,2,0.5,2))
# plot.peak(df.peak = hela_NAT,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
plot.peak(df.peak = hela_NAT.fix,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

# par(mar=c(0.5,2,0.5,2))
# plot.peak(df.peak = hela_g_hic_TAD.overlap.raw.fix,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF0000",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)

par(mar=c(0.5,2,0.5,2))
plot.TAD(n_hic_matrix.part,ylim=c(-100,0),maxBound = 0.958,mat.upper = F)

par(mar=c(0.5,2,0.5,2))
plot.barplot(hela_CTCF.fix,chrom_index,chrom_start,chrom_end)


quantile(n_hic_matrix.part,prob=0.958)
quantile(g_hic_matrix.part,prob=0.97)


dev.off()

par(mar=c(1,1,1,1))

rm(list=ls())


##############################################################################
# Fig3 make NAT other region
##############################################################################
chrom_end = NULL
for(i in c(1:23)){
  chrom_len = chrom.length(i)
  chrom_end = c(chrom_end,chrom_len-1)
}

## make hg19 chromsome bed format file
hg19_genome_table = data.frame(chrom_name = c(paste("chr",c(1:22),sep = ""),"chrX"),
                               chrom_start = rep(0,23),
                               chrom_end = chrom_end,
                               name = rep("chromsome",23),
                               value = c(1:23))

## NAT complement
genome_table = NULL

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  
  
  chrom_table = hg19_genome_table[hg19_genome_table$chrom_name==chrom_name,]
  chrom_NAT_table = hela_NAT[hela_NAT$V1==chrom_name,]
  
  start_vector = c(0,chrom_NAT_table$V3 + 1)
  end_vector = c(chrom_NAT_table[1,2] - 1,chrom_NAT_table$V2[2:length(chrom_NAT_table$V2)]-1,chrom.length(chrom_index)-1)
  
  chrom_region_table = data.frame(chrom_name = rep(chrom_name,length(start_vector)),
                                  start = start_vector,
                                  end = end_vector,
                                  info = sprintf("region_%d",c(1:length(start_vector))),
                                  value = c(1:length(start_vector))
                                  )
  genome_table = rbind(genome_table,chrom_region_table)
}

write.table(genome_table,file = "~/menghw_HD/data_table/Hela_NAT_complement.bed",col.names = F,row.names = F,quote = F,sep = "\t")


## NAD complement
genome_table = NULL

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  
  
  chrom_table = hg19_genome_table[hg19_genome_table$chrom_name==chrom_name,]
  chrom_NAD_table = hela_NAD[hela_NAD$V1==chrom_name,]
  
  start_vector = c(0,chrom_NAD_table$V3 + 1)
  end_vector = c(chrom_NAD_table[1,2] - 1,chrom_NAD_table$V2[2:length(chrom_NAD_table$V2)]-1,chrom.length(chrom_index)-1)
  
  chrom_region_table = data.frame(chrom_name = rep(chrom_name,length(start_vector)),
                                  start = start_vector,
                                  end = end_vector,
                                  info = sprintf("region_%d",c(1:length(start_vector))),
                                  value = c(1:length(start_vector))
  )
  genome_table = rbind(genome_table,chrom_region_table)
}

genome_table.fix = genome_table[c(-79,-80),]
genome_table.fix.fix = genome_table.fix[c(-20),]

write.table(genome_table.fix.fix,file = "~/menghw_HD/data_table/Hela_NAD_complement.bed",col.names = F,row.names = F,quote = F,sep = "\t")


##############################################################################
# Fig3 plot chromosome and insulation score
##############################################################################
load(file="~/menghw_HD/R_code/my_function/MyPlot_v05.RData")

hela_NAT = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
hela_NAD = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep = "\t")
hela_n_TAD_score = read.table("~/menghw_HD/data_table/TAD_region/hela-n-hic-merged.insulation_clean_up.bed",header = F,sep = "\t")
hela_n_TAD_score$V4[is.na(hela_n_TAD_score$V4)]  = 0

chrom_index = 8
binsize = 10e3
image.path = sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_matrix_chr_%d_%d.RData",chrom_index,binsize)
load(file=image.path)

n_hic_matrix_ActD.raw = read.table("~/menghw_HD/our_data/Hela-nucleolus-hic-actD-hg19/matrix/raw_matrix/MAPQ20/binsize_10000/chr_8_10000_MAPQ20.txt",header = F,sep = ",")

g_hic_matrix.z.score.value = sum(g_hic_matrix.norm[upper.tri(g_hic_matrix.norm,diag = F)]) / 10e6
n_hic_matrix.z.score.value = sum(n_hic_matrix[upper.tri(n_hic_matrix,diag = F)]) / 10e6
n_hic_matrix_ActD.z.score.value = sum(n_hic_matrix_ActD.raw[upper.tri(n_hic_matrix_ActD.raw,diag = F)]) / 10e6

chrom_start = 128.5e6
chrom_end = 128.8e6

g_hic_matrix.part = matrix.part(g_hic_matrix.norm,chrom_index,chrom_start,chrom_end,binsize) / g_hic_matrix.z.score.value
n_hic_matrix.part = matrix.part(n_hic_matrix,chrom_index,chrom_start,chrom_end,binsize) / n_hic_matrix.z.score.value
n_hic_matrix_ActD.part = matrix.part(n_hic_matrix_ActD.raw,chrom_index,chrom_start,chrom_end,binsize) / n_hic_matrix_ActD.z.score.value
differ_matrix.part = n_hic_matrix.part - g_hic_matrix.part

quantile(g_hic_matrix.part,prob=0.95)
quantile(n_hic_matrix.part,prob=0.95)
quantile(n_hic_matrix_ActD.part,prob=0.95)




# png(file="~/menghw_HD/test7.png",width = 1000,height = 3000)

fig.facet <- layout(matrix(c(1:8),nrow = 8,byrow = FALSE),width = rep(10,4),heights = c(5,3,5,3,5,3,3))
layout.show(fig.facet)

par(mar=c(0.5,4,0.5,1))
plot.TAD(g_hic_matrix.part,mat.part = F,maxBound = 0.4)
axis(side=1)

par(mar=c(0.5,4,0.5,1))
plot(colSums(g_hic_matrix.part),type="l",frame.plot =F,xaxt="n",xlab="",ylab="",lwd=2.5,ylim=c(0,2.5e4),yaxt="n",col="#0F414F")
axis(side=2,at=c(0,2.5e4),lwd=3,col="#022A35")
abline(h=mean(colSums(g_hic_matrix.part)),lwd=1,pch=10,lty=2,col="black")

par(mar=c(0.5,4,0.5,1))
plot.TAD(n_hic_matrix.part,mat.part = F,ylim = c(-50,50),maxBound = 0.95)
axis(side=1)

par(mar=c(0.5,4,0.5,1))
plot(colSums(n_hic_matrix.part),type="l",frame.plot =F,xaxt="n",xlab="",ylab="",lwd=2.5,ylim=c(0,4e4),col="#0F414F",yaxt="n")
axis(side=2,at=c(0,4e4),lwd=3,col="#022A35")
abline(h=mean(colSums(n_hic_matrix.part)),lwd=1,pch=10,lty=2,col="black")

par(mar=c(0.5,4,0.5,1))
plot.TAD(n_hic_matrix_ActD.part,mat.part = F,ylim = c(-50,50),maxBound = 0.95)

par(mar=c(0.5,4,0.5,1))
chrom_name = chrom.name(chrom_index)
chrom_TAD_score = hela_n_TAD_score[hela_n_TAD_score$V1==chrom_name,]
chrom_TAD_score.x = chrom_TAD_score[chrom_TAD_score$V2>=chrom_start & chrom_TAD_score$V2<=chrom_end,2]
chrom_TAD_score.vector = chrom_TAD_score[chrom_TAD_score$V2>=chrom_start & chrom_TAD_score$V2<=chrom_end,4]
plot(x=chrom_TAD_score.x,y=chrom_TAD_score.vector,type="l",frame.plot =F,xaxt="n",xlab="",ylab="",lwd=2.5,ylim=c(-3,3),col="#0F414F",yaxt="n")
axis(side=2,at=c(-3,3),lwd=3,col="#022A35")
abline(h=0,lwd=1,pch=10,lty=2,col="black")
# par(mar=c(0,2,0.5,2))
# plot(colSums(n_hic_matrix.part),type="l",frame.plot =F,xaxt="n",xlab="",ylab="",lwd=2.5,ylim=c(0,100000))

# par(mar=c(0.5,4,0.5,1))
# plot(colSums(differ_matrix.part),type="l",frame.plot =F,xaxt="n",xlab="",ylab="",lwd=2.5,ylim=c(-2.5e4,2.5e4))

dev.off()


hela_all_table = read.table("~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep = "\t")
hela_all_table.fix = hela_all_table[hela_all_table$g_count>0,]

NAD_vector = hela_all_table.fix$n_count / hela_all_table.fix$g_count
hela_NAD = hela_all_table.fix[NAD_vector>2,]

hela_all_table.fix = hela_all_table[hela_all_table$g_hic>0,]
NAT_vector = hela_all_table.fix$n_hic / hela_all_table.fix$g_hic
hela_NAT = hela_all_table.fix[NAT_vector >2,]



hela_NAT[hela_NAT$chr_index=="chr16",]
