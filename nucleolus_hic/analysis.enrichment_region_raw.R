# Howard MENG on 2016-July-1st
# make enrichment region statistics
# We deal peak calling with MACS2 software combine with --broad parameter

# # 函数定义区域
# chrom.name <- function(chrom_index = 1){
#   if(chrom_index <=22){
#     chrom_name = paste0("chr",chrom_index)
#   }else if(chrom_index == 23){
#     chrom_name = "chrX"
#   }else if(chrom_index == 24){
#     chrom_name = "chrY"
#   }else if (chrom_index == 25){
#     chrom_name = "rDNA"
#   }
#   return(chrom_name)
# }
# 
# # hg19 长度
HG19_LEN.total = 3095720411
HG19_LEN.n = 234350292
HG19_LEN.not_n = 2861370119

# # A/B compartment table 
# binsize = 1e5
# hela_ab_table = read.table("~/menghw_HD/data_table/Hela-genome-hic_AB_all_100000_ice.table",header = T,sep = "\t")
# 
# # Hela-nucleolus-hic enrichment table
# hela_enrich_table_real = read.table("~/menghw_HD/data_table/Hela-nucleolus-hic_enrich_BroadPeak.bed",header = F,sep = "\t")
# hela_enrich_table_real = hela_enrich_table_real[,c(1:5)]
# colnames(hela_enrich_table_real) = c("chrom_name","start","end","peak_name","value")
# 
# ## hic enrichment region windowsize=100Kb
# hela_enrich_table_combine = read.table("~/menghw_HD/data_table/Hela-n-hic_enrich_filter_w1e5.table",header = T)
# 
# # N-seq and G-seq call peak
# ## with MACS2 --broad parameter and merge the peak within a window size of 100Kb
# ## 读取NAD区域
hela_NAD_table = read.table("~/menghw_HD/data_table/raw_table/Hela-n-seq_enrich_BroadPeak_control.bed",header = F,sep = "\t")
hela_NAD_table.fix.vector = hela_NAD_table$V5 >= 50 
hela_NAD_table.fix = hela_NAD_table[hela_NAD_table.fix.vector,]


# save.image("~/menghw_HD/R_image/analysis_enrichment.RData")

# 载入数据, 数据内容由上面的部分生成
rm(list=ls())
load(file = "~/menghw_HD/R_image/analysis_enrichment.RData")

binsize = 1e5

# A/B compartment length
HG19_LEN.compartment = nrow(hela_ab_table) * binsize 

# A/B ration in whole genome
a_ratio = sum(hela_ab_table$value[hela_ab_table$info=="A"]) / nrow(hela_ab_table)
b_ratio = sum(hela_ab_table$value[hela_ab_table$info=="B"]) / nrow(hela_ab_table)

# E/NE ration in whole genome
e_combine.ratio = sum(hela_enrich_table_combine$end - hela_enrich_table_combine$start) / HG19_LEN.total
e_real.ratio = sum(hela_enrich_table_real$end - hela_enrich_table_real$start) / HG19_LEN.total

# hela-hic-enrich.fix
quantile(hela_enrich_table_real$value,prob=seq(0,1,0.05))
# hela_enrich_table_real.fix = hela_enrich_table_real[hela_enrich_table_real$value>=120,]
hela_enrich_table_real.fix = read.table("~/menghw_HD/data_table/Hela-n-hic_enrich_filter_80_w1e4.table",header = F,sep = "\t")
colnames(hela_enrich_table_real.fix) = colnames(hela_enrich_table_real)
e_ratio = sum(hela_enrich_table_real.fix$end - hela_enrich_table_real.fix$start) / HG19_LEN.total

# NAD region in whole genome
NAD.ratio = sum(hela_NAD_table$end - hela_NAD_table$start) / HG19_LEN.total

# 计算overlap
## overlap between A/B and E/NE
genome_b_ratio.vector = rep(0,23)
for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_ab_table = hela_ab_table[hela_ab_table$chrom_name==chrom_name,]
  genome_b_ratio.vector[chrom_index] = nrow(chrom_ab_table[chrom_ab_table$info=="A",]) / nrow(chrom_ab_table)
}

e_b_ratio.vector = rep(0,23)
for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_ab_table = hela_ab_table[hela_ab_table$chrom_name==chrom_name,]
  # chrom_enrich_table = hela_enrich_table_combine[hela_enrich_table_combine$chrom_name==chrom_name,]
  # chrom_enrich_table = hela_enrich_table_real[hela_enrich_table_real$chrom_name==chrom_name,]
  chrom_enrich_table = hela_enrich_table_real.fix[hela_enrich_table_real.fix$chrom_name==chrom_name,]
  
  chrom_overlab_table = data.frame(matrix(0,nrow = nrow(chrom_enrich_table),3))
  colnames(chrom_overlab_table) = c("A_length","B_length","Total_length")
  
  for(i in c(1:nrow(chrom_enrich_table))){
    peak.start = chrom_enrich_table$start[i]
    peak.end = chrom_enrich_table$end[i]
    peak.start.index = peak.start %/% binsize
    peak.end.index = peak.end %/% binsize
    
    for(index in c(peak.start.index:peak.end.index)){
      region.start = index * binsize
      region.end = region.start + binsize
      region.info = chrom_ab_table$info[chrom_ab_table$start==region.start]
      if(length(region.info) >0){
        add.length.list = sort(c(region.start,region.end,peak.start,peak.end))
        add.length = add.length.list[3] - add.length.list[2]
        if(region.info == "A"){
          chrom_overlab_table$A_length[i] = chrom_overlab_table$A_length[i] + add.length
        }else{
          chrom_overlab_table$B_length[i] = chrom_overlab_table$B_length[i] + add.length
        }
        add.length = 0
      }  
    }
  }
  e_b_ratio.vector[chrom_index] = colSums(chrom_overlab_table)[2] / sum(chrom_overlab_table[,c(1:2)])
}

# 计算非enrichment region
ne_b_ratio.vector = rep(0,23)
for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_ab_table = hela_ab_table[hela_ab_table$chrom_name==chrom_name,]
  chrom_enrich_table = hela_enrich_table_real.fix[hela_enrich_table_real.fix$chrom_name==chrom_name,]
  
  chrom_overlab_table = data.frame(matrix(0,nrow = nrow(chrom_enrich_table),3))
  colnames(chrom_overlab_table) = c("A_length","B_length","Total_length")
  
  ne_region_index.vector = c(0)
  
  e_region_index.start = chrom_enrich_table$start %/% binsize
  e_region_index.end = chrom_enrich_table$end %/% binsize
  e_region_index = c(unique(c(e_region_index.start,e_region_index.end)))
  
  ne_region_filter = chrom_ab_table$start != e_region_index[1]*binsize
  
  for(index in e_region_index){
    ne_region_filter = ne_region_filter &  (chrom_ab_table$start != index*binsize)
  }
  
  chrom_ne_region_table = chrom_ab_table[ne_region_filter,]
  ne_b_ratio.vector[chrom_index] = nrow(chrom_ne_region_table[chrom_ne_region_table$info=="B",])/nrow(chrom_ne_region_table)
}

# B-ratio boxplot
g_color = rgb(181,134,84,maxColorValue = 255)
e_color = rgb(0,87,55,maxColorValue = 255)
boxplot(ne_b_ratio.vector,genome_b_ratio.vector,e_b_ratio.vector,col=c("#6495ED","#FF8C00","#FF4500"),lwd=3,cex.axis=2)
box(lwd=3)

# plot RPM
hela_table = read.table("~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep="\t")
filter_vector = hela_table$not_n_rate>=0.5 & hela_table$copy_number!=0
hela_table.fix = hela_table[filter_vector,]

## 计算各染色体的RPM
e_RPM.vector = rep(0,23)
ne_RPM.vector = rep(0,23)

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_table = hela_table.fix[hela_table.fix$chr_index==chrom_name,]
  chrom_enrich_table = hela_enrich_table_real.fix[hela_enrich_table_real.fix$chrom_name==chrom_name,]
  
  
  ne_region_index.vector = c(0)
  
  e_region_index.start = chrom_enrich_table$start %/% binsize
  e_region_index.end = chrom_enrich_table$end %/% binsize
  e_region_index = c(unique(c(e_region_index.start,e_region_index.end)))
  
  ne_region_filter = chrom_table$region_min != e_region_index[1]*binsize
  
  for(index in e_region_index){
    ne_region_filter = ne_region_filter &  (chrom_table$region_min != index*binsize)
  }
  
  e_RPM.vector[chrom_index] = sum((chrom_table$n_count / chrom_table$copy_number)[e_region_index] ,na.rm = T) / (length(e_region_index) - length(which(is.na(chrom_table$copy_number[e_region_index])))) / binsize * 1e6 
  ne_RPM.vector[chrom_index] = sum((chrom_table$n_count / chrom_table$copy_number)[ne_region_filter] ,na.rm = T) / (length(ne_region_filter) - length(which(is.na(chrom_table$copy_number[ne_region_filter])))) / binsize * 1e6 
}

g_color = rgb(181,134,84,maxColorValue = 255)
e_color = rgb(0,87,55,maxColorValue = 255)
barplot(c(ne_RPM.vector,e_RPM.vector),col =c("#6495ED","#FF8C00"))
boxplot(ne_RPM.vector,e_RPM.vector,col=c("#6495ED","#FF8C00","#FF4500"),lwd=3,cex.axis=2)
box(lwd=3)


chrom_index = 1
chrom_name = chrom.name(chrom_index)
RKM = hela_table$


# plot Enrichment region
non_e_ratio = 1 - NAD.ratio - e_ratio
pie(c(non_e_ratio,NAD.ratio,e_ratio),col = c("#00BFFF","#FF8C00","#FF6347"),border = F,labels = c("Not-NAIR","NAD","NAIR"))
pie(c(non_e_ratio,NAD.ratio,e_ratio),col = c("#00BFFF","#FF8C00","#FF6347"),border = F,labels = "")


# plot NAD / G-Band / NAIR
load(file="~/menghw_HD/R_image/MyPlot.RData")

hela_g_band = read.table(file = "~/menghw_HD/reference/hg19_G-band.txt",header = F,sep = "\t")
colnames(hela_g_band) = c("chrom_name","start","end","name","value","info")

chrom_index = 1
chrom_start = 0
chrom_end = chrom.length(chrom_index)


png(file = "~/menghw_HD/R_image/Hela-hic-figure/test.png",width = 500,height = 5000)

fig.facet <- layout(matrix(c(1:91),nrow = 91,byrow = FALSE),width = rep(10,91),heights = c(rep(c(1,1,1,0.5),22),c(1,1,1)))
layout.show(fig.facet)

for(chrom_index in (1:23)){
  chrom_start = 0
  chrom_end = chrom.length(1)
  # plot NAD
  par(mar=c(0.5,1,0,1))
  plot.peak(df.peak = hela_NAD_table,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#303030",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  # plot NAIR
  par(mar=c(0.5,1,0,1))
  plot.peak(df.peak = hela_enrich_table_real.fix,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#EE5C42",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  # plot G-Band
  par(mar=c(0.5,1,0,1))
  plot.peak(hela_g_band,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom.length(chrom_index),color_pattern = T,track_color = "black",track_pattern = T,xaxt = F,yaxt = F,rect_border = T,rect_border.lwd = 2,xylwd = 3,xcex = 3)
  
  par(mar=c(0,1,0,1))
  #plot(x=c(1),y=c(1),type="n",frame.plot = F,ylab = "",xlab = "",xaxt = "n",yaxt="n")
}

dev.off()




  
# write.table(hela_enrich_table_real.fix,file = "~/menghw_HD/data_table/Hela-n-hic_enrich_filter_80.table",quote = F,col.names = T,row.names = F,sep = "\t")


for(i in c(1:nrow(chrom_enrich_table))){
  peak.start = chrom_enrich_table$start[i]
  peak.end = chrom_enrich_table$end[i]
  peak.start.index = peak.start %/% binsize
  peak.end.index = peak.end %/% binsize
  
  for(index in c(peak.start.index:peak.end.index)){
    region.start = index * binsize
    region.end = region.start + binsize
    region.info = chrom_ab_table$info[chrom_ab_table$start==region.start]
    if(length(region.info) >0){
      add.length.list = sort(c(region.start,region.end,peak.start,peak.end))
      add.length = add.length.list[3] - add.length.list[2]
      if(region.info == "A"){
        chrom_overlab_table$A_length[i] = chrom_overlab_table$A_length[i] + add.length
      }else{
        chrom_overlab_table$B_length[i] = chrom_overlab_table$B_length[i] + add.length
      }
      add.length = 0
    }  
  }
}



chrom_name = chrom.name(chrom_index)
chrom_ab_table = hela_ab_table[hela_ab_table$chrom_name==chrom_name,]
chrom_enrich_table = hela_enrich_table_real.fix[hela_enrich_table_real.fix$chrom_name==chrom_name,]

chrom_overlab_table = data.frame(matrix(0,nrow = nrow(chrom_enrich_table),3))
colnames(chrom_overlab_table) = c("A_length","B_length","Total_length")

ne_region_index.vector = c(0)

e_region_index.start = chrom_enrich_table$start %/% binsize
e_region_index.end = chrom_enrich_table$end %/% binsize
e_region_index = c(unique(c(e_region_index.start,e_region_index.end)))

ne_region_filter = chrom_ab_table$start != e_region_index[1]*binsize

for(index in e_region_index){
  ne_region_filter = ne_region_filter &  (chrom_ab_table$start != index*binsize)
}

chrom_ne_region_table = chrom_ab_table[ne_region_filter,]
ne_b_ratio.vector[chrom_index] = nrow(chrom_ne_region_table[chrom_ne_region_table$info=="B",])/nrow(chrom_ne_region_table)


# a function for counting overlap bewteen two bed files
overlap <- function(table_1,table_2,level_vetor=c("chr1")){
  
  overlap_table = data.frame(row.names = c("chrom_name","start","end","name","start.1","end.1","name.1","start.2","end.2","name.2"))
  
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
                                name.2 = region_2.name
                                )
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

hela_enrich_table_real.fix = hela_enrich_table_real[hela_enrich_table_real$value>=500,]

x = overlap(hela_ab_table,hela_enrich_table_real.fix)

View(x)


load("~/menghw_HD/R_image/MyPlot.RData")

# Hela NAD 区域
hela_n_seq_control = read.table("~/menghw_HD/data_table/raw_table/Hela-n-seq_enrich_BroadPeak_control.bed",sep = "\t")
hela_n_seq_no_control = read.table("~/menghw_HD/data_table/raw_table/Hela-n-seq_no_control.broadPeak",sep="\t")
## 总长度
hela_n_seq_control.len = sum(hela_n_seq_control$V3 - hela_n_seq_control$V2)
hela_n_seq_no_control.len = sum(hela_n_seq_no_control$V3 - hela_n_seq_no_control$V2)
## overlap
## figure


# Hela hic enrichment区域 
hela_n_hic_control.cis = read.table("~/menghw_HD/data_table/raw_table/Hela-n-hic_control_cis.broadPeak",sep = "\t")
hela_n_hic_control.all = read.table("~/menghw_HD/data_table/raw_table/Hela-n-hic_enrich_BroadPeak_control_all.bed",sep = "\t")

hela_n_hic_no_control.cis = read.table("~/menghw_HD/data_table/raw_table/Hela-n-hic_no_control_cis.broadPeak",sep="\t")
hela_n_hic_no_control.all = read.table("~/menghw_HD/data_table/raw_table/Hela-n-hic_no_control_all.broadPeak",sep="\t")

# fix
# quantile(hela_n_hic_control.cis$V5,prob=seq(0,1,0.05))
# hela_n_hic_control.cis.fix = hela_n_hic_control.cis[hela_n_hic_control.cis$V5>=100,]

quantile(hela_n_hic_control.all$V5,prob=seq(0,1,0.05))
hela_n_hic_control.all.fix = hela_n_hic_control.all[hela_n_hic_control.all$V5>=60,]

quantile(hela_n_hic_no_control.all$V5,prob=seq(0,1,0.05))
quantile(hela_n_hic_no_control.all$V3 - hela_n_hic_no_control.all$V2,prob=seq(0,1,0.05))
hela_n_hic_no_control.all.fix = hela_n_hic_no_control.all[hela_n_hic_no_control.all$V5>=60,]
hela_n_hic_no_control.all.fix = hela_n_hic_no_control.all[hela_n_hic_no_control.all$V5>=60 & ((hela_n_hic_no_control.all$V3-hela_n_hic_no_control.all$V2)>=3000),]


# quantile(hela_n_hic_no_control.cis$V5,prob=seq(0,1,0.05))
# hela_n_hic_no_control.cis.fix = hela_n_hic_no_control.cis[hela_n_hic_no_control.cis$V5>=115,]

## 总长度
hela_n_hic_control.all.len = sum(hela_n_hic_control.all$V3 - hela_n_hic_control.all$V2)
hela_n_hic_no_control.all.len = sum(hela_n_hic_no_control.all$V3 - hela_n_hic_no_control.all$V2)

hela_n_hic_control.all.fix.len = sum(hela_n_hic_control.all.fix$V3 - hela_n_hic_control.all.fix$V2)
hela_n_hic_no_control.all.fix.len = sum(hela_n_hic_no_control.all.fix$V3 - hela_n_hic_no_control.all.fix$V2)

## 长度分布
hist(hela_n_hic_control.all$V3 - hela_n_hic_control.all$V2)
hist(hela_n_hic_no_control.all$V3 - hela_n_hic_no_control.all$V2)
ks.test(x=(hela_n_hic_control.all$V3 - hela_n_hic_control.all$V2),y=(hela_n_hic_no_control.all$V3 - hela_n_hic_no_control.all$V2))
## overlap


## figure
png("~/menghw_HD/R_image/test.png",width = 800,height = 4000)
fig.nrow = 46
fig.ncol = 1
fig.facet <- layout(matrix(c(1:fig.nrow),nrow = fig.nrow,byrow = FALSE),width = rep(10,fig.nrow),heights = rep(1,46))
layout.show(fig.facet)


par(mar=c(1,1,1,1))
plot.peak(hela_n_hic_control.all,chrom_index = chrom_index,chrom_len = chrom.length(chrom_index),track_pattern = T,color_pattern = F,xaxt = F,yaxt = F)
plot.peak(hela_n_hic_no_control.all,chrom_index = chrom_index,chrom_len = chrom.length(chrom_index),track_pattern = T,color_pattern = F,xaxt = T,yaxt = F)

for(chrom_index in c(1:23)){
  par(mar=c(0.5,1,0.5,1))
  plot.peak(hela_n_hic_seq_control,chrom_index = chrom_index,chrom_len = chrom.length(chrom_index),track_pattern = T,color_pattern = F,xaxt = F,yaxt = F)
  # plot.peak(hela_n_hic_no_control.all.fix,chrom_index = chrom_index,chrom_len = chrom.length(chrom_index),track_pattern = T,color_pattern = F,xaxt = F,yaxt = F)
  par(mar=c(0.5,1,0,1))
  plot.peak(hela_LAD_region.fix,chrom_index = chrom_index,chrom_len = chrom.length(chrom_index),track_pattern = T,color_pattern = F,xaxt = F,yaxt = F)
  # plot.peak(hela_n_hic_control.all.fix,chrom_index = chrom_index,chrom_len = chrom.length(chrom_index),track_pattern = T,color_pattern = F,xaxt = F,yaxt = F)  
}


plot.peak(hela_n_hic_no_control.all.fix,chrom_index = chrom_index,chrom_len = chrom.length(chrom_index),track_pattern = T,color_pattern = F,xaxt = T,yaxt = F)

plot.peak(hg19_gene_table,chrom_index = chrom_index,chrom_len = chrom.length(chrom_index),track_pattern = T,color_pattern = F,xaxt = T,yaxt = F)
dev.off()

hg19_gene_table = read.table("~/menghw_HD/reference/gene_gtf/hg19_annotation.gene.bed",header = T,sep = "\t")
View(hg19_gene_table[1:10,])

hela_n_hic_seq_control = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-hic/seq-hic-control/hela-n-hic_peaks.broadPeak",sep = "\t",header = F)
hela_n_hic_seq_control_broad = read.table("~/menghw_HD/our_data/Data_PeakCalling/with_control/hela-n-hic/seq-hic-control/hela-n-hic-broad_peaks.broadPeak",sep = "\t",header = F)

hela_n_hic_seq_control_broad = read.table("~/menghw_HD/data_table/raw_table/Hela-n-hic_seq_control_w1e4.bed",sep = "\t",header = F)


quantile(hela_n_hic_seq_control$V3 - hela_n_hic_seq_control$V2,prob=seq(0,1,0.05))
quantile(hela_n_hic_seq_control$V5,prob=seq(0,1,0.05))

hela_n_seq.filtervector = ((hela_n_hic_seq_control$V3 - hela_n_hic_seq_control$V2) >= 1000) & (hela_n_hic_seq_control$V5 >= 30)

hela_n_hic_seq_control.fix = hela_n_hic_seq_control[hela_n_seq.filtervector,]


sum(hela_n_hic_seq_control.fix$V3 - hela_n_hic_seq_control.fix$V2)
write.table(hela_n_hic_seq_control.fix,file = "~/menghw_HD/data_table/fix_table/Hela-n-hic_seq_control_fix.table",quote = F,col.names = F,row.names = F,sep = "\t")

hela_n_hic_enrich = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_seq_control_fix_w1e4.bed",header = F,sep = "\t")
sum(hela_n_hic_enrich$V3 - hela_n_hic_enrich$V2)
quantile(hela_n_hic_enrich$V3 - hela_n_hic_enrich$V2)


hist(hela_n_hic_enrich$V3 - hela_n_hic_enrich$V2,breaks = seq(1e3,9e6,by=2e4),xlim = c(0,5e5))

par(mar=c(5,5,2,2))
hist(log10(hela_n_hic_enrich$V3 - hela_n_hic_enrich$V2),breaks = seq(3,7,0.1),xlab = "Log10 Length",ylab = "Frequency",main="",cex.axis=2,cex.lab=2,lwd=3,col="#FF8247",)


hist(hela_n_hic_seq_control$V3 - hela_n_hic_seq_control$V2,breaks = seq(0,3e5,by=5e2),xlim = c(0,2e4))


# enrichment region length 
par(mar=c(5,5,2,2))
hist(log10(hela_n_hic_enrich$V3 - hela_n_hic_enrich$V2),breaks = seq(3,7,0.1),xlab = "Log10 Length",ylab = "Frequency",main="",cex.axis=2,cex.lab=2,lwd=3,col="#FF8247",)
abline(v=median(log10(hela_n_hic_enrich$V3 - hela_n_hic_enrich$V2)),lwd=5,pch=10,lty=2,col="black")
text(x=4,y=800,"median=3.98",pos=4,cex=2.5)

# enrichment ration among chrom
load(file="~/menghw_HD/R_image/MyPlot.RData")
enrichment_rate.vector = c(1:23)
for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  enrichment_rate.vector[chrom_index] = sum((hela_n_hic_enrich$V3 - hela_n_hic_enrich$V2)[hela_n_hic_enrich$V1==chrom_name]) / chrom.length(chrom_index)
}


chrom_name.vector = c(1:23)
for(chrom_index in c(1:23)){
  chrom_name.vector[chrom_index] = chrom.name(chrom_index)
}

nair_rate = sum(hela_n_hic_enrich$V3 - hela_n_hic_enrich$V2) / HG19_LEN.total

par(mar=c(5,5,5,2))
barplot(enrichment_rate.vector*100,ylim = c(0,75),border = NA,space =0.5 ,col = ifelse(enrichment_rate.vector>=nair_rate,"#FF8247","#7EC0EE"),xaxt="n",yaxt="n")
# par(new=T)
axis(side=2,at=c(0,25,round(nair_rate*100,1),50,75),labels = paste(c(0,25,round(nair_rate*100,1),50,75),"%",sep = "") ,cex.axis=3,lwd = 5)
text(x = seq(0.3,34,1.5),y= -1,srt=-75, pos=4, xpd=TRUE,labels=chrom_name.vector,cex=2)
abline(h=nair_rate*100,lty=2,lwd=5,col="black")

# plot enrichment region 
fig.facet <- layout(matrix(c(1:46),nrow = 46,byrow = FALSE),width = rep(10,46),heights = rep(1,46))
layout.show(fig.facet)

for(chrom_index in (1:23)){
  chrom_start = 0
  chrom_end = chrom.length(1)
  # plot NAD
  par(mar=c(0.5,1,0,1))
  plot.peak(df.peak = hela_NAD_table.fix,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#303030",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  # plot NAIR
  par(new=T)
  plot.peak(df.peak = hela_n_hic_enrich,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF8247",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)
  
  # plot G-Band
  par(mar=c(0.5,1,0,1))
  plot.peak(hg19_g_band_track,chrom_index,chrom_start = chrom_start,chrom_end = chrom_end,chrom_len = chrom.length(chrom_index),color_pattern = T,track_color = "black",track_pattern = T,xaxt = F,yaxt = F,rect_border = T,rect_border.lwd = 2,xylwd = 3,xcex = 3)
  
  # par(mar=c(0,1,0,1))
  #plot(x=c(1),y=c(1),type="n",frame.plot = F,ylab = "",xlab = "",xaxt = "n",yaxt="n")
}
plot.peak(hg19_g_band_track,chrom_index = 1,chrom_len = chrom.length(chrom_index),chrom_start = 0,chrom_end = chrom.length(chrom_index),color_pattern = T,track_pattern = T,track_color = "black",xaxt = F,yaxt = F,rect_border = T)
par(new=T)
plot.peak(hela_n_hic_enrich,chrom_index = 1,chrom_len = chrom.length(chrom_index),chrom_start = 0,chrom_end = chrom.length(chrom_index),color_pattern = T,track_pattern = T,track_color = "black",xaxt = F,yaxt = F,rect_border = T)



