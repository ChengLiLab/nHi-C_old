rm(list=ls())

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
    chrom_index = as.numeric(sub("chr","",chrom_name))
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
load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")
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


####################################
## matrix normalization
####################################
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
hela_all_table = read.table("~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep = "\t")
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
relation_matrix[4,2] = relation_matrix[2,4]

relation_matrix[3,4] = cor(x=hela_all_table$n_count,y=hela_all_table$n_hic)
relation_matrix[4,3] = relation_matrix[3,4]

colnames(relation_matrix) = c("WGS","Hi-C","NS","nHi-C")
rownames(relation_matrix) = colnames(relation_matrix)

par(mar=c(1,1,1,1))
plot.matrix(relation_matrix-0.1,max_bound = 1)

library(pheatmap)
library(gplots)
par(mar=c(2,2,2,2))
pheatmap(relation_matrix,color = colorpanel(100,low = "yellow",high = "red"),cluster_rows = F,cluster_cols = F,cex=1.5,display_numbers=T,fontsize_number = 18,number_color = "black")






