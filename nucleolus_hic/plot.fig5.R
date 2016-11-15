#############################################################################################################
# 2016-09-10 Howard MENG
# plot after ActD treatment 
#############################################################################################################

################################################################
# plot all-all matrix
################################################################
# load all-all matrix data 
# save.image(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-all-all-matrix_1000000.RData")
rm(list=ls())

load(file="~/menghw_HD/R_image/Hela-hic-Rimage/Hela-all-all-matrix_1000000.RData")
load(file="~/menghw_HD/R_code/my_function/MyPlot_v07.RData")
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")

g_hic_all_matrix.fix = g_hic_all_matrix.norm / hic_count.hela_g_hic.total * 300e6
quantile(g_hic_all_matrix.fix,prob=c(0,0.97))
png(filename = "~/menghw_HD/R_image/Hela-hic-figure/Hela-g-hic_all_all_matrix_v2.png",width = 3000,height = 3000)
par(mar=c(1,1,1,1))
# plot.matrix(g_hic_all_matrix.fix,bound.max = 0.99)
plot.matrix(log2(g_hic_all_matrix.fix + 1),bound.max = 0.97)
dev.off()


n_hic_all_matrix.fix = n_hic_all_matrix / hic_count.hela_n_hic.total * 300e6
quantile(n_hic_all_matrix.fix,prob=c(0,0.97))
png(filename = "~/menghw_HD/R_image/Hela-hic-figure/Hela-n-hic_all_all_matrix_v2.png",width = 3000,height = 3000)
par(mar=c(1,1,1,1))
# plot.matrix(n_hic_all_matrix.fix,bound.max = 0.99)
plot.matrix(log2(n_hic_all_matrix.fix + 1),bound.max = 0.97)
dev.off()


g_hic_ActD_all_matrix.fix = g_hic_ActD_all_matrix.norm / hic_count.hela_g_hic_ActD.total * 300e6
quantile(g_hic_ActD_all_matrix.fix,prob=c(0,0.97))
png(filename = "~/menghw_HD/R_image/Hela-hic-figure/Hela-g-hic-ActD_all_all_matrix_v2.png",width = 3000,height = 3000)
par(mar=c(1,1,1,1))
# plot.matrix(g_hic_ActD_all_matrix.fix,bound.max = 0.99)
plot.matrix(log2(g_hic_ActD_all_matrix.fix + 1),bound.max = 0.97)
dev.off()


n_hic_ActD_all_matrix.fix = n_hic_ActD_all_matrix / hic_count.hela_n_hic_ActD.total * 300e6
quantile(n_hic_ActD_all_matrix.fix,prob=c(0,0.97))
png(filename = "~/menghw_HD/R_image/Hela-hic-figure/Hela-n-hic-ActD_all_all_matrix_v2.png",width = 3000,height = 3000)
par(mar=c(1,1,1,1))
# plot.matrix(n_hic_ActD_all_matrix.fix,bound.max = 0.99)
plot.matrix(log2(n_hic_ActD_all_matrix.fix + 1),bound.max = 0.97)
dev.off()


################################################################
# differ matrix 
################################################################
# g-hic-ActD / g-hic

g_hic_all_matrix.fix = g_hic_all_matrix / hic_count.hela_g_hic.total * 300e6
g_hic_ActD_all_matrix.fix = g_hic_ActD_all_matrix / hic_count.hela_g_hic_ActD.total * 300e6
g_hic_differ_matrix = g_hic_ActD_all_matrix.fix / g_hic_all_matrix.fix

g_hic_differ_matrix.log = log2(g_hic_differ_matrix)
g_hic_differ_matrix.log[g_hic_all_matrix.fix==0 | g_hic_differ_matrix==0] = 0 

quantile(g_hic_differ_matrix.log,prob=c(0.0005,0.98))
png(filename = "~/menghw_HD/R_image/Hela-hic-figure/Hela-g-hic_differ_all_all_matrix_log_v2.png",width = 3000,height = 3000)
plot.matrix(g_hic_differ_matrix.log,col.min = "blue",col.max = "red",bound.max = 0.98,bound.min = 0.0005,col.boundary = 0,n_block_color = "#888888")
dev.off()

# n-hic-ActD / n-hic
n_hic_differ_matrix = n_hic_ActD_all_matrix.fix / n_hic_all_matrix.fix

n_hic_differ_matrix.log = log2(n_hic_differ_matrix)
n_hic_differ_matrix.log[n_hic_all_matrix.fix==0 | n_hic_differ_matrix==0] = 0 

quantile(n_hic_differ_matrix.log,prob=c(0.00005,1))

png(filename = "~/menghw_HD/R_image/Hela-hic-figure/Hela-n-hic_differ_all_all_matrix_log_v4.png",width = 3000,height = 3000)
plot.matrix(n_hic_differ_matrix.log,col.min = "blue",col.max = "red",bound.max = 1,bound.min = 0.0005,col.boundary = 0,n_block_color = "#888888")
dev.off()


################################################################
# plot matrix in one chromosome
################################################################
rm(list=ls())
chrom_index = 19
binsize = 50e3

load(file=sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_all_matrix_chr_%d_%d.RData",chrom_index,binsize))
load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")
hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
colnames(hela_NAIR_table) = c("chrom_name","start","end","peak_index","info")

chrom_name = chrom.name(chrom_index)
# conver as.matrix
segment.row = c(30e6,59e6)

# segment.row = c(0,chrom.length(chrom_index))
segment.col = segment.row
chrom_start = segment.col[1]
chrom_end = segment.col[2]

g_hic_matrix.part = matrix.part.corner(g_hic_matrix,chrom_index,segment.row,segment.col,binsize)
n_hic_matrix.part = matrix.part.corner(n_hic_matrix,chrom_index,segment.row,segment.col,binsize)
g_hic_ActD_matrix.part = matrix.part.corner(g_hic_ActD_matrix,chrom_index,segment.row,segment.col,binsize)
n_hic_ActD_matrix.part = matrix.part.corner(n_hic_ActD_matrix,chrom_index,segment.row,segment.col,binsize)


# make matrix with the same sequencing depth
g_hic_matrix.fix = g_hic_matrix.part / hic_count.hela_g_hic.total * 300e6
n_hic_matrix.fix = n_hic_matrix.part / hic_count.hela_n_hic.total * 300e6
g_hic_ActD_matrix.fix = g_hic_ActD_matrix.part / hic_count.hela_g_hic_ActD.total * 300e6
n_hic_ActD_matrix.fix = n_hic_ActD_matrix.part / hic_count.hela_n_hic_ActD.total * 300e6

######### make differ matrix #########
# g/g and n/n 
g_hic_differ_matrix =  g_hic_ActD_matrix.fix / g_hic_matrix.fix
g_hic_differ_matrix.log = log2(g_hic_differ_matrix)
g_hic_differ_matrix.log[g_hic_matrix.fix==0 | g_hic_differ_matrix==0] = 0
quantile(g_hic_differ_matrix.log,prob=c(0.035,0.98))

# hela-g-hic-ActD / hela-g-hic
fig.facet <- layout(matrix(c(1:2),nrow = 2,byrow = FALSE),width = c(10,10),heights = c(10,1))
layout.show(fig.facet)
par(mar=c(0,1,1,1))
plot.matrix(g_hic_differ_matrix.log,col.min = "blue",col.max = "red",bound.max = 0.98,bound.min = 0.035,col.boundary = 0,n_block_color = "white")
par(mar=c(0.5,1,0,1))
plot.peak(df.peak = hela_NAIR_table,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)


# hela-n-hic-ActD / hela-n-hic
n_hic_differ_matrix =  n_hic_ActD_matrix.fix / n_hic_matrix.fix
n_hic_differ_matrix.log = log2(n_hic_differ_matrix)
n_hic_differ_matrix.log[n_hic_matrix.fix==0 | n_hic_differ_matrix==0] = 0
quantile(n_hic_differ_matrix.log,prob=c(0.027,0.95))

fig.facet <- layout(matrix(c(1:2),nrow = 2,byrow = FALSE),width = c(10,10),heights = c(10,1))
layout.show(fig.facet)

par(mar=c(0,1,1,1))
plot.matrix(n_hic_differ_matrix.log,col.min = "blue",col.max = "red",bound.max = 0.95,bound.min = 0.027,col.boundary = 0,n_block_color = "white")

par(mar=c(0.5,1,0,1))
plot.peak(df.peak = hela_NAIR_table,chrom_index,chrom_len,chrom_start,chrom_end,yaxt = F,xaxt = F,axis_x_len=axis_x_len,track_color = "#FF7304",track_pattern = T,color_pattern = F,peak.lwd = 1.5,xylwd = 3)


# n/g and n/g ActD 
g_hic_differ_matrix.log = log2(g_hic_differ_matrix)
g_hic_differ_matrix.log[g_hic_all_matrix.fix==0 | g_hic_differ_matrix==0] = 0 

quantile(g_hic_differ_matrix.log,prob=c(0.0005,0.98))
png(filename = "~/menghw_HD/R_image/Hela-hic-figure/Hela-g-hic_differ_all_all_matrix_log_v2.png",width = 3000,height = 3000)
plot.matrix(g_hic_differ_matrix.log,col.min = "blue",col.max = "red",bound.max = 0.98,bound.min = 0.0005,col.boundary = 0,n_block_color = "#888888")


# plot diff matrix in one chromosome 



################################################################
# fig4 cis/trans ratio
################################################################
rm(list=ls())
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")
hic_count_table = data.frame(hela_g_hic = c(hic_count.hela_g_hic.cis,hic_count.hela_g_hic.trans),
                             hela_g_hic_ActD = c(hic_count.hela_g_hic_ActD.cis,hic_count.hela_g_hic_ActD.trans),
                             hela_n_hic = c(hic_count.hela_n_hic.cis,hic_count.hela_n_hic.trans),
                             hela_n_hic_ActD = c(hic_count.hela_n_hic_ActD.cis,hic_count.hela_n_hic_ActD.trans))
rownames(hic_count_table) = c("Cis_count","Trans_count")

barplot.matrix = rbind(as.matrix(hic_count_table)[1,] / colSums(as.matrix(hic_count_table)),as.matrix(hic_count_table)[2,] / colSums(as.matrix(hic_count_table)))

par(mar=c(8,6,2,1),family="Arial",font=1)
par(mar=c(2,4,2,2),family="Arial",font=1)
barplot(barplot.matrix,col = c("#C936D3","#89A110"),border = F,axes = F,axisnames=F)
axis(side=2,at=seq(0,1,0.25),labels = F,cex.axis=3,lwd=3,tck=-0.1)
text(x=-0.9,y=seq(0,1,0.25),labels = paste0(seq(0,100,25),"%"),cex = 2,xpd=T,pos = NULL)
text(x=c(1.2,2.2,3.4,4.6),y=-0.05,labels = c("Hi-C","Hi-C ActD","nHi-C","nHi-C ActD"),cex = 2,xpd=T,srt = 45,pos = 2)



################################################################
# A/B changes
################################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v07.RData")
hela_ab_value_table = read.table("~/menghw_HD/data_table/Hela_AB_compart_value_fix.table",header = T,sep = "\t")
# ActD_ice_coef = c(0,0,0,0,0,0,-1,1,1,-1,-1,-1,0,0,-1,-1,-1,-1,-1,1,-1,1,0)
# hela_ab_value_table = cbind(hela_ab_value_table,ActD_ice_coef)
# write.table(hela_ab_value_table,file = "~/menghw_HD/data_table/Hela_AB_compart_value_fix.table",col.names = T,row.names = F,quote = F,sep = "\t")
chrom_index_list = c(7:12,15:22)

# 1. 有多少A/B 发生变化了
hela_g_hic_ActD_AB_table = read.table("~/menghw_HD/data_table/Hela-g-hic-ActD_AB.table",header = T,sep = "\t")
hela_g_hic_ActD_AB_table$start = hela_g_hic_ActD_AB_table$start - 100e3
hela_g_hic_ActD_AB_table$end = hela_g_hic_ActD_AB_table$end - 100e3

hela_g_hic_AB_table = read.table("~/menghw_HD/data_table/Hela-genome-hic_AB_all_100000_ice.table",header = T,sep = "\t")
hela_g_hic_AB_table.part = hela_g_hic_AB_table[hela_g_hic_AB_table$chrom_name %in% chrom.name.factor(chrom_index_list),]

# A-B & B-A & B-B & A-A
AB_same_table = NULL
AB_diff_table = NULL

for(chrom_index in chrom_index_list){
  chrom_same_table = NULL
  chrom_diff_table = NULL
  chrom_name = chrom.name(chrom_index)
  
  print(chrom_name)
  chrom_AB_table.ctrl = hela_g_hic_AB_table.part[hela_g_hic_AB_table.part$chrom_name==chrom_name,]
  chrom_AB_table.ActD = hela_g_hic_ActD_AB_table[hela_g_hic_ActD_AB_table$chrom_name==chrom_name,]
  
  for(i in c(1:nrow(chrom_AB_table.ActD))){
    ActD_start = as.integer(chrom_AB_table.ActD$start[i])
    ActD_type = as.character(chrom_AB_table.ActD$type[i])
    ctrl_type = as.character(chrom_AB_table.ctrl$info[chrom_AB_table.ctrl$start==ActD_start])
    
    if(length(ctrl_type)==0){
      chrom_diff_table = rbind(chrom_diff_table,chrom_AB_table.ActD[i,])
    }else if(ActD_type == ctrl_type){
      chrom_same_table = rbind(chrom_same_table,chrom_AB_table.ActD[i,])
    }else{
      chrom_diff_table = rbind(chrom_diff_table,chrom_AB_table.ActD[i,])
    }
  }
  
  AB_same_table = rbind(AB_same_table,chrom_same_table)
  AB_diff_table = rbind(AB_diff_table,chrom_diff_table)
}

write.table(AB_diff_table,file = "~/menghw_HD/data_table/Hela-g-hic-ActD_AB_diff.table",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(AB_same_table,file = "~/menghw_HD/data_table/Hela-g-hic-ActD_AB_same.table",col.names = T,row.names = F,quote = F,sep = "\t")

# 2. A/B 变化
A2A_rate = length(which(AB_same_table$type == "A")) / nrow(hela_g_hic_ActD_AB_table)
B2B_rate = length(which(AB_same_table$type == "B")) / nrow(hela_g_hic_ActD_AB_table)
A2B_rate = length(which(AB_diff_table$type == "B")) / nrow(hela_g_hic_ActD_AB_table)
B2A_rate = length(which(AB_diff_table$type == "A")) / nrow(hela_g_hic_ActD_AB_table)

par(mar=c(1,1,1,1))
pie.data = c(A2A_rate,B2B_rate,A2B_rate,B2A_rate)
# pie(pie.data,col = c("#000000AA","#00000088","#4573D5","#FFBE40"),cex=2,border = F,labels = c("A to A","B to B","A to B","B to A"))
pie(pie.data,col = c("#000000AA","#00000088","#06276F","#FFBE40"),cex=2,border = F,labels = "")

# 2. A/B 变化能解释多少变异的基因
load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")
load(file="~/menghw_HD/R_code/my_function/overlap_gene.RData")
hela_gene_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/Hela_gene_exp_all_rmdup.table",header = T,sep = "\t")
hela_gene_exp_table = hela_gene_exp_table[hela_gene_exp_table$gene_type=="protein_coding",]
chrom_index_list = c(7:12,15:22)
chrom_name_list = chrom.name.factor(chrom_index_list)

hela_gene_exp_table.part = NULL
for(chrom_index in chrom_index_list){
  chrom_name = chrom.name(chrom_index)
  print(chrom_name)
  chrom_gene_exp_table = hela_gene_exp_table[hela_gene_exp_table$chrom_name==chrom_name,]
  hela_gene_exp_table.part = rbind(hela_gene_exp_table.part, chrom_gene_exp_table)
}

hela_gene_exp_table.part.up = hela_gene_exp_table.part[hela_gene_exp_table.part$fold_change>0 & hela_gene_exp_table.part$significant=="yes",]
hela_gene_exp_table.part.down = hela_gene_exp_table.part[hela_gene_exp_table.part$fold_change<0 & hela_gene_exp_table.part$significant=="yes",]
hela_gene_exp_table.part.equal = hela_gene_exp_table.part[hela_gene_exp_table.part$test_state=="OK" & hela_gene_exp_table.part$significant=="no",]

A2B_table = AB_diff_table[AB_diff_table$type == "B",]
B2A_table = AB_diff_table[AB_diff_table$type == "A",]

A2A_table = AB_same_table[AB_same_table$type == "A",]
B2B_table = AB_same_table[AB_same_table$type == "B",]


A2B_gene_up_overlap = overlap.gene(A2B_table,hela_gene_exp_table.part,chrom_name_list)
B2A_gene_up_overlap = overlap.gene(B2A_table,hela_gene_exp_table.part,chrom_name_list)

A2A_gene_up_overlap = overlap.gene(A2A_table,hela_gene_exp_table.part,chrom_name_list)
B2B_gene_up_overlap = overlap.gene(B2B_table,hela_gene_exp_table.part,chrom_name_list)

save.image(file="~/menghw_HD/data_table/R_data/AB_gene_exp_overlap.RData")
# 3. boxplot
A2A_gene_up_overlap.test = A2A_gene_up_overlap[A2A_gene_up_overlap$test_state=="OK",]
B2A_gene_up_overlap.test = B2A_gene_up_overlap[B2A_gene_up_overlap$test_state=="OK" & B2A_gene_up_overlap$significant=="yes",]
A2B_gene_up_overlap.test = A2B_gene_up_overlap[A2B_gene_up_overlap$test_state=="OK" & A2B_gene_up_overlap$significant=="yes",]
B2B_gene_up_overlap.test = B2B_gene_up_overlap[B2B_gene_up_overlap$test_state=="OK",]


# p value < 0.01 
A2B_gene_up_overlap.test = A2B_gene_up_overlap[A2B_gene_up_overlap$test_state=="OK" & A2B_gene_up_overlap$significant=="yes" & A2B_gene_up_overlap$p_value <= 0.01,]
A2B_gene_up_overlap.test = unique(A2B_gene_up_overlap.test)

B2B_gene_up_overlap.test = B2B_gene_up_overlap[B2B_gene_up_overlap$test_state=="OK" & B2B_gene_up_overlap$significant=="yes" & B2B_gene_up_overlap$p_value <= 0.01,]
B2B_gene_up_overlap.test = unique(B2B_gene_up_overlap.test)

A2A_vector = as.vector(A2A_gene_up_overlap.test$fold_change)
B2B_vector = as.vector(B2B_gene_up_overlap.test$fold_change) 
A2B_vector = as.vector(A2B_gene_up_overlap.test$fold_change) 
B2A_vector = as.vector(B2A_gene_up_overlap.test$fold_change)

par(mar=c(3,3,1,1),family="Arial",font=1)
# boxplot(c(A2A_vector,B2B_vector),A2B_vector,B2A_vector,ylim=c(-6,6),col=c("#00000088","#06276F","#FFBE40"),xaxt="n",yaxt="n",frame.plot=F)
boxplot(c(A2A_vector,B2B_vector),A2B_vector,B2A_vector,ylim=c(-6,6),col=c("#00000088","#06276F","#FFBE40"),xaxt="n",yaxt="n",frame.plot=F)
axis(side=1,at=c(1:3),labels =F,outer = F,cex.axis=2.5,lwd=3,tck=-0.05)
# text(x=c(1:3),y=-7.5,labels = c("No change","A to B","B to A"),cex = 2,xpd=T)

axis(side=2,labels = F,cex.axis=3,lwd=3,tck=-0.05,at=seq(-6,6,3))
text(x=0.1,y=seq(-6,6,3),labels = seq(-6,6,3),cex = 2,xpd=T)



##################################################################################################
# interaction-distance plot
##################################################################################################
rm(list=ls())

######################################################
# function
######################################################
# matrix to table 
load(file="~/menghw_HD/R_code/my_function/MyPlot_v07.RData")

make.dist.table <- function(hic_matrix,binsize){
  # 输入chromosome hic matrix 与binsize 生成distance-interaction dataframe
  mat.nrow = dim(hic_matrix)[1]
  mat.ncol = dim(hic_matrix)[2]
  
  count_vector = rep(0,mat.nrow)
  for(row.index in c(1:mat.nrow)){
    x = c(row.index : mat.nrow)
    y = c(1:length(x))
    inter_count = 0
    for(xy.index in c(1:length(x))){
      x0 = x[xy.index]
      y0 = y[xy.index]
      inter_count = inter_count + hic_matrix[x0,y0]
    }
    count_vector[row.index] = inter_count
  }
  
  # 生成返回的data frame
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }else if(chrom_index == 25){
    chrom_name = "rDNA"
  }
  
  start_vector = as.integer(c(1:mat.nrow) * binsize)
  
  inter_dist_table = data.frame(chrom_name = rep(chrom_name,length(count_vector)),
                                start = start_vector,
                                hic_count = count_vector)
  return(inter_dist_table)
}


######################################################
# plot interaction-distance 
######################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/make_dist_table.RData")
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")


binsize = 100e3

for(chrom_index in c(1:23)){
  print(chrom_index)
  
  load(file=sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_all_matrix_chr_%d_%d.RData",chrom_index,binsize))
  
  # conver as.matrix
  g_hic_matrix.raw = as.matrix(g_hic_matrix)
  n_hic_matrix.raw = as.matrix(n_hic_matrix)
  g_hic_ActD_matrix.raw = as.matrix(g_hic_ActD_matrix)
  n_hic_ActD_matrix.raw = as.matrix(n_hic_ActD_matrix)
  
  # make matrix with the same sequencing depth
  g_hic_matrix.fix = g_hic_matrix.raw / hic_count.hela_g_hic.total * 300e6
  n_hic_matrix.fix = n_hic_matrix.raw / hic_count.hela_n_hic.total * 300e6
  g_hic_ActD_matrix.fix = g_hic_ActD_matrix.raw / hic_count.hela_g_hic_ActD.total * 300e6
  n_hic_ActD_matrix.fix = n_hic_ActD_matrix.raw / hic_count.hela_n_hic_ActD.total * 300e6
  
  # make inter-dist table
  g_hic_dist_table = make.dist.table(g_hic_matrix.fix,binsize)
  g_hic_ActD_dist_table = make.dist.table(g_hic_ActD_matrix.fix,binsize)
  n_hic_dist_table = make.dist.table(n_hic_matrix.fix,binsize)
  n_hic_ActD_dist_table = make.dist.table(n_hic_ActD_matrix.fix,binsize)
  
  ylim = c(0,7)
  xlim = c(5,8.5)
  
  png(file=sprintf("~/menghw_HD/R_image/Hela-hic-figure/inter_dist_plot/hic_inter_dist_chr_%d_%d_v3.png",chrom_index,binsize),width = 1000,height = 1000)
  
  plot(x=log10(g_hic_dist_table$start),y=log10(g_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#1533AD",frame.plot=F,yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  axis(side=2,at=c(0:7))
  
  par(new=T)
  plot(x=log10(g_hic_ActD_dist_table$start),y=log10(g_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#6F81D6",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  par(new=T)
  plot(x=log10(n_hic_dist_table$start),y=log10(n_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FF8B00",frame.plot=F,yaxt="n",xlab="",xaxt="n",ylab="",xlim=xlim,ylim=ylim)
  
  par(new=T)
  plot(x=log10(n_hic_ActD_dist_table$start),y=log10(n_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FFC640",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  dev.off()
  
  
}



######################################################
# NAIR 区的region plot inter-dist  
######################################################
rm(list=ls())
hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
colnames(hela_NAIR_table) = c("chrom_name","start","end","peak_index","info")

# fix function 
chrom_index = 18
binsize = 100e3
load(file=sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_all_matrix_chr_%d_%d.RData",chrom_index,binsize))

hic_matrix = as.matrix(g_hic_matrix)
NAIR_table = hela_NAIR_table


make.NAIR.dist.table <- function(hic_matrix,binsize,NAIR_table,chrom_name){
  # 输入chromosome hic matrix 与binsize 生成distance-interaction dataframe
  chrom_NAIR_table = NAIR_table[NAIR_table[,1]==chrom_name,]
  colnames(chrom_NAIR_table) = c("chrom_name","start","end","peak_index","info")
  NAIR_bin_index_list = NULL
  for(i in c(1:nrow(chrom_NAIR_table))){
    bin.index.start = chrom_NAIR_table$start[i] %/% binsize 
    bin.index.end = chrom_NAIR_table$end[i] %/% binsize
    NAIR_bin_index_list = c(NAIR_bin_index_list,c(bin.index.start:bin.index.end))
  }
  
  mat.nrow = dim(hic_matrix)[1]
  mat.ncol = dim(hic_matrix)[2]
  
  count_vector = rep(0,mat.nrow)
  for(row.index in c(1:mat.nrow)){
    x = c(row.index : mat.nrow)
    y = c(1:length(x))
    inter_count = 0
    for(xy.index in c(1:length(x))){
      x0 = x[xy.index]
      y0 = y[xy.index]
      if(x0 %in% NAIR_bin_index_list){
        inter_count = inter_count + hic_matrix[x0,y0]  
      }
    }
    count_vector[row.index] = inter_count
  }
  
  # 生成返回的data frame
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }else if(chrom_index == 25){
    chrom_name = "rDNA"
  }
  
  start_vector = as.integer(c(1:mat.nrow) * binsize)
  
  inter_dist_table = data.frame(chrom_name = rep(chrom_name,length(count_vector)),
                                start = start_vector,
                                hic_count = count_vector)
  return(inter_dist_table)
}


###### plot 
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v07.RData")
load(file="~/menghw_HD/R_code/my_function/make_dist_table.RData")
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")
hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
colnames(hela_NAIR_table) = c("chrom_name","start","end","peak_index","info")
binsize = 100e3

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  print(chrom_name)
  
  load(file=sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_all_matrix_chr_%d_%d.RData",chrom_index,binsize))
  
  # conver as.matrix
  g_hic_matrix.raw = as.matrix(g_hic_matrix)
  n_hic_matrix.raw = as.matrix(n_hic_matrix)
  g_hic_ActD_matrix.raw = as.matrix(g_hic_ActD_matrix)
  n_hic_ActD_matrix.raw = as.matrix(n_hic_ActD_matrix)
  
  # make matrix with the same sequencing depth
  g_hic_matrix.fix = g_hic_matrix.raw / hic_count.hela_g_hic.total * 300e6
  n_hic_matrix.fix = n_hic_matrix.raw / hic_count.hela_n_hic.total * 300e6
  g_hic_ActD_matrix.fix = g_hic_ActD_matrix.raw / hic_count.hela_g_hic_ActD.total * 300e6
  n_hic_ActD_matrix.fix = n_hic_ActD_matrix.raw / hic_count.hela_n_hic_ActD.total * 300e6
  
  # make inter-dist table
  g_hic_dist_table = make.NAIR.dist.table(g_hic_matrix.fix,binsize,hela_NAIR_table,chrom_name)
  g_hic_ActD_dist_table = make.NAIR.dist.table(g_hic_ActD_matrix.fix,binsize,hela_NAIR_table,chrom_name)
  n_hic_dist_table = make.NAIR.dist.table(n_hic_matrix.fix,binsize,hela_NAIR_table,chrom_name)
  n_hic_ActD_dist_table = make.NAIR.dist.table(n_hic_ActD_matrix.fix,binsize,hela_NAIR_table,chrom_name)
  
  ylim = c(0,7)
  xlim = c(0,chrom.length(chrom_index))
  
  png(file=sprintf("~/menghw_HD/R_image/Hela-hic-figure/inter_dist_plot/chr_%d_%d_NAIR_dist.png",chrom_index,binsize),width = 1000,height = 1000)
  
  plot(x=g_hic_dist_table$start,y=log10(g_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#1533AD",frame.plot=F,yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  axis(side=2,at=c(0:7))
  
  par(new=T)
  plot(x=g_hic_ActD_dist_table$start,y=log10(g_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#6F81D6",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  par(new=T)
  plot(x=n_hic_dist_table$start,y=log10(n_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FF8B00",frame.plot=F,yaxt="n",xlab="",xaxt="n",ylab="",xlim=xlim,ylim=ylim)
  
  par(new=T)
  plot(x=n_hic_ActD_dist_table$start,y=log10(n_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FFC640",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  dev.off()
  
}


#################################################################
# NAIR 区的region plot inter-dist  但是使用的是已经算好的RData
#################################################################
rm(list=ls())

binsize = 100e3

for(chrom_index in c(1:23)){
  print(chrom_index)
  load(sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela_inter_dist/hela_inter_dist_chr_%d_%d.RData",chrom_index,binsize))
  
  ylim = c(0,7)
  xlim = c(5,8.5)
  
  png(file=sprintf("~/menghw_HD/R_image/Hela-hic-figure/inter_dist_plot/hic_inter_dist_NAIR_chr_%d_%d_v3.png",chrom_index,binsize),width = 1000,height = 1000)
  
  plot(x=log10(g_hic_dist_table$start),y=log10(g_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#1533AD",frame.plot=F,yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  axis(side=2,at=c(0:7))
  
  par(new=T)
  plot(x=log10(g_hic_ActD_dist_table$start),y=log10(g_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#6F81D6",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  par(new=T)
  plot(x=log10(n_hic_dist_table$start),y=log10(n_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FF8B00",frame.plot=F,yaxt="n",xlab="",xaxt="n",ylab="",xlim=xlim,ylim=ylim)
  
  par(new=T)
  plot(x=log10(n_hic_ActD_dist_table$start),y=log10(n_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=3,col="#FFC640",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)
  
  title(main = sprintf("chr%d NAIR",chrom_index))
  
  dev.off()
  
}


######################################################
# NAIR-NAIR interaction v.s. NAIR-self interaction
######################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/make_dist_table.RData")
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")



######################################################
# plot chr19 inter-dist plot 精修格式
######################################################
rm(list=ls())

chrom_index = 19
binsize = 100e3

load(file="~/menghw_HD/R_code/my_function/make_dist_table.RData")
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")
load(file=sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_all_matrix_chr_%d_%d.RData",chrom_index,binsize))
load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")

# segment.row = c(25e6,95e6)
segment.row = c(0,chrom.length(chrom_index))
segment.col = segment.row
chrom_start = segment.col[1]
chrom_end = segment.col[2]

g_hic_matrix.raw = matrix.part.corner(g_hic_matrix,chrom_index,segment.row,segment.col,binsize)
n_hic_matrix.raw = matrix.part.corner(n_hic_matrix,chrom_index,segment.row,segment.col,binsize)
g_hic_ActD_matrix.raw = matrix.part.corner(g_hic_ActD_matrix,chrom_index,segment.row,segment.col,binsize)
n_hic_ActD_matrix.raw = matrix.part.corner(n_hic_ActD_matrix,chrom_index,segment.row,segment.col,binsize)

# make matrix with the same sequencing depth
g_hic_matrix.fix = g_hic_matrix.raw / hic_count.hela_g_hic.total * 300e6
n_hic_matrix.fix = n_hic_matrix.raw / hic_count.hela_n_hic.total * 300e6
g_hic_ActD_matrix.fix = g_hic_ActD_matrix.raw / hic_count.hela_g_hic_ActD.total * 300e6
n_hic_ActD_matrix.fix = n_hic_ActD_matrix.raw / hic_count.hela_n_hic_ActD.total * 300e6

# make inter-dist table
g_hic_dist_table = make.dist.table(g_hic_matrix.fix,binsize)
g_hic_ActD_dist_table = make.dist.table(g_hic_ActD_matrix.fix,binsize)
n_hic_dist_table = make.dist.table(n_hic_matrix.fix,binsize)
n_hic_ActD_dist_table = make.dist.table(n_hic_ActD_matrix.fix,binsize)

ylim = c(0,7)
xlim = c(5,8)

# png(file=sprintf("~/menghw_HD/R_image/Hela-hic-figure/inter_dist_plot/chr_inter_dist_%d_%d_v3.png",chrom_index,binsize),width = 1000,height = 1000)

line.lwd = 8
# g-hic 
plot(x=log10(g_hic_dist_table$start),y=log10(g_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=line.lwd,col="#1533AD",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)

# g-hic-ActD
par(new=T)
plot(x=log10(g_hic_ActD_dist_table$start),y=log10(g_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=line.lwd,col="#6F81D6",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)

axis(side=2,at=c(0:7),labels = F,lwd=3,tck=-0.05)
axis(side=1,at=seq(5,8,0.5),labels = F,lwd=3,tck=-0.05)


# par(mar=c(2,2,0,0))
# n-hic 
plot(x=log10(n_hic_dist_table$start),y=log10(n_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=line.lwd,col="#FF8B00",frame.plot=F,yaxt="n",xlab="",xaxt="n",ylab="",xlim=xlim,ylim=ylim)

# n-hic-ActD
par(new=T)
plot(x=log10(n_hic_ActD_dist_table$start),y=log10(n_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=line.lwd,col="#FFC640",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)

axis(side=2,at=c(0:7),labels = F,lwd=3,tck=-0.05)
axis(side=1,at=seq(5,8,0.5),labels = F,lwd=3,tck=-0.05)



######################################################
# plot chr19 inter-dist plot 精修格式 NAIR
######################################################
rm(list=ls())

chrom_index = 19
binsize = 100e3
load(sprintf("~/menghw_HD/R_image/Hela-hic-Rimage/Hela_inter_dist/hela_inter_dist_chr_%d_%d.RData",chrom_index,binsize))

ylim = c(0,7)
xlim = c(5,8)

# png(file=sprintf("~/menghw_HD/R_image/Hela-hic-figure/inter_dist_plot/hic_inter_dist_NAIR_chr_%d_%d_v3.png",chrom_index,binsize),width = 1000,height = 1000)
line.lwd = 8

plot(x=log10(g_hic_dist_table$start),y=log10(g_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=line.lwd,col="#1533AD",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)

par(new=T)
plot(x=log10(g_hic_ActD_dist_table$start),y=log10(g_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=line.lwd,col="#6F81D6",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)

par(new=T)
plot(x=log10(n_hic_dist_table$start),y=log10(n_hic_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=line.lwd,col="#FF8B00",frame.plot=F,yaxt="n",xlab="",xaxt="n",ylab="",xlim=xlim,ylim=ylim)

par(new=T)
plot(x=log10(n_hic_ActD_dist_table$start),y=log10(n_hic_ActD_dist_table$hic_count),type="l",pch=19,cex=0.3,lwd=line.lwd,col="#FFC640",frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=xlim,ylim=ylim)

axis(side=2,at=c(0:7),labels = F,lwd=3,tck=-0.05)
axis(side=1,at=seq(5,8,0.5),labels = F,lwd=3,tck=-0.05)

dev.off()



