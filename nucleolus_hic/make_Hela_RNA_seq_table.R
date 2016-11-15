##################################################################################################
# 整理Hela-RNA-Seq data 
##################################################################################################

###############################################
# loading data
###############################################
# RNA_seq_table.ActD.1 = read.table("~/menghw_HD/data_table/Hela-RNA-seq/cufflink_data/Hela-RNA-seq_gene_FPKM_ActD_rep1.table",header = T,sep = '\t')
# RNA_seq_table.ActD.2 = read.table("~/menghw_HD/data_table/Hela-RNA-seq/cufflink_data/Hela-RNA-seq_gene_FPKM_ActD_rep2.table",header = T,sep = '\t')
# RNA_seq_table.ActD.3 = read.table("~/menghw_HD/data_table/Hela-RNA-seq/cufflink_data/Hela-RNA-seq_gene_FPKM_ActD_rep3.table",header = T,sep = '\t')
# 
# RNA_seq_table.ctrl.1 = read.table("~/menghw_HD/data_table/Hela-RNA-seq/cufflink_data/Hela-RNA-seq_gene_FPKM_ctrl_rep1.table",header = T,sep = '\t')
# RNA_seq_table.ctrl.2 = read.table("~/menghw_HD/data_table/Hela-RNA-seq/cufflink_data/Hela-RNA-seq_gene_FPKM_ctrl_rep2.table",header = T,sep = '\t')
# RNA_seq_table.ctrl.3 = read.table("~/menghw_HD/data_table/Hela-RNA-seq/cufflink_data/Hela-RNA-seq_gene_FPKM_ctrl_rep3.table",header = T,sep = '\t')
# 
# RNA_seq_table.total = read.table("~/menghw_HD/data_table/Hela-RNA-seq/raw/Hela_exp_fixed.table",header = T,sep = '\t')
# 
# RNA_seq_table.merge = RNA_seq_table.total[,c(1:11,15,16)]
# RNA_seq_table.merge = cbind(RNA_seq_table.merge,rep(0,nrow(RNA_seq_table.merge)),rep(0,nrow(RNA_seq_table.merge)),rep(0,nrow(RNA_seq_table.merge)),rep(0,nrow(RNA_seq_table.merge)),rep(0,nrow(RNA_seq_table.merge)),rep(0,nrow(RNA_seq_table.merge)))
# colnames(RNA_seq_table.merge) = c(colnames(RNA_seq_table.merge)[1:(ncol(RNA_seq_table.merge)-6)],"ctrl_1","ctrl_2","ctrl_3","ActD_1","ActD_2","ActD_3")
# 
# nrow(RNA_seq_table.total[RNA_seq_table.total$gene_type=="protein_coding",])
# ctrl_1_table = RNA_seq_table.ctrl.1[-grep(pattern = "ctrl.",RNA_seq_table.ctrl.1$gene_id),]
# 
# 
# 
# RNA_seq_table.raw = read.table("~/menghw_HD/data_table/Hela-RNA-seq/Hela_gene_exp.diff",header = T,sep = "\t")
# RNA_seq_table.raw.fix = RNA_seq_table.raw[RNA_seq_table.raw$gene!="-",]
# 
# hg19_GTF_table = read.table("~/menghw_HD/reference/gene_gtf/hg19_coding.gene.bed",header = F,sep = '\t')



###############################################
# make table
###############################################
# for(row.index in c(1:nrow(RNA_seq_table.merge))){
#   if(row.index %% 1000 == 0){
#     print(row.index / nrow(RNA_seq_table.merge) * 100)
#   }
#   
#   gene_name = as.character(RNA_seq_table.merge$gene_name[row.index])
#   
#   ## ctrl
#   ctrl_1 = RNA_seq_table.ctrl.1$FPKM[grep(pattern = gene_name,RNA_seq_table.ctrl.1$gene_short_name)]
#   if(length(ctrl_1) > 1){
#     print(ctrl_1)
#   }
#   # ctrl_1 = ctrl_1[1]
#   # 
#   # ctrl_2 = RNA_seq_table.ctrl.2$FPKM[grep(pattern = gene_name,RNA_seq_table.ctrl.2$gene_short_name)]
#   # ctrl_2 = ctrl_2[1]
#   # 
#   # ctrl_3 = RNA_seq_table.ctrl.3$FPKM[grep(pattern = gene_name,RNA_seq_table.ctrl.3$gene_short_name)]
#   # ctrl_3 = ctrl_3[1]
#   # 
#   # ## ActD
#   # actd_1 = RNA_seq_table.ActD.1$FPKM[grep(pattern = gene_name,RNA_seq_table.ActD.1$gene_short_name)]
#   # actd_1 = actd_1[1]
#   # 
#   # actd_2 = RNA_seq_table.ActD.2$FPKM[grep(pattern = gene_name,RNA_seq_table.ActD.1$gene_short_name)]
#   # actd_2 = actd_2[1]
#   # 
#   # actd_3 = RNA_seq_table.ActD.3$FPKM[grep(pattern = gene_name,RNA_seq_table.ActD.1$gene_short_name)]
#   # actd_3 = actd_3[1]
#   # 
#   # RNA_seq_table.merge$ctrl_1[row.index] = ctrl_1
#   # RNA_seq_table.merge$ctrl_2[row.index] = ctrl_2
#   # RNA_seq_table.merge$ctrl_3[row.index] = ctrl_3
#   # RNA_seq_table.merge$ActD_1[row.index] = actd_1
#   # RNA_seq_table.merge$ActD_2[row.index] = actd_2
#   # RNA_seq_table.merge$ActD_3[row.index] = actd_3
# }
# 
# 
# write.table(RNA_seq_table.merge,file="~/menghw_HD/data_table/Hela-RNA-seq/Hela-RNA-seq.table",col.names = T,row.names = F,quote = F,sep = "\t")
# 
# ###############################################
# # plot heatmap
# ###############################################
# RNA_seq_table.merge.fix = RNA_seq_table.merge[RNA_seq_table.merge$test_state=="OK" ,c(4,14:19)]
# 
# RNA_seq_table.merge.fix = na.omit(RNA_seq_table.merge.fix)
# 
# RNA_seq_table.merge.fix.mat = as.matrix(RNA_seq_table.merge.fix[,c(2:7)])
# rownames(RNA_seq_table.merge.fix.mat) = RNA_seq_table.merge.fix$gene_name
# 
# library(pheatmap)
# RNA_seq_table.merge.fix.mat.log = log2(RNA_seq_table.merge.fix.mat)
# RNA_seq_table.merge.fix.mat.log[RNA_seq_table.merge.fix.mat==0] = 0
# 
# pheatmap(mat =log2(RNA_seq_table.merge.fix.mat+1),border_color = F,legend_breaks=c(0:5))
# 
# dev.off()




###############################################
# 对test state == OK 的数据进行补齐
###############################################
hela_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/Hela-RNA-seq_fixed_v2.table",header = T,sep = "\t")

hela_exp_table.fix = hela_exp_table[hela_exp_table$gene_type=="protein_coding" & (!is.na(hela_exp_table$gene_type)),]

hela_exp_table.fix = hela_exp_table[! is.na(hela_exp_table$gene_type),]
hela_exp_table.fix = hela_exp_table.fix[hela_exp_table.fix$test_state=="OK",]


for(row.index in c(1:nrow(hela_exp_table.fix))){
  for(col.index in c(14:16)){
    if(is.na(hela_exp_table.fix[row.index,col.index])){
      if(hela_exp_table.fix$FPKM.control[row.index]==0){
        hela_exp_table.fix[row.index,col.index] = 0
      }else{
        hela_exp_table.fix[row.index,col.index] = abs(rnorm(n = 1,mean = hela_exp_table.fix$FPKM.control[row.index] ,sd = sqrt(hela_exp_table.fix$FPKM.control[row.index])))  
      }
      
    }
  }
  
  for(col.index in c(17:19)){
    if(is.na(hela_exp_table.fix[row.index,col.index])){
      if(hela_exp_table.fix$FPKM.ActD[row.index]==0){
        hela_exp_table.fix[row.index,col.index] = 0  
      }else{
        hela_exp_table.fix[row.index,col.index] = abs(rnorm(n = 1,mean =hela_exp_table.fix$FPKM.ActD[row.index] ,sd = sqrt(hela_exp_table.fix$FPKM.ActD[row.index])))  
      }
    }
  }
}

write.table(hela_exp_table.fix,file = "~/menghw_HD/data_table/Hela-RNA-seq/Hela_gene_exp_all.table",col.names = T,row.names = F,sep = "\t",quote = F)



###############################################
# plot heatmap 判断数据重复性
###############################################

library(pheatmap)
library(gplots)

RNA_seq_matrix = as.matrix(hela_exp_table.fix[,c(14:19)])
colnames(RNA_seq_matrix) = c(paste("ctrl_",c(1:3),sep = ""),paste("ActD_",c(1:3),sep = ""))
rownames(RNA_seq_matrix) = hela_exp_table.fix$gene_name

RNA_seq_matrix[RNA_seq_matrix==0] = 1
RNA_seq_matrix.log = log2(RNA_seq_matrix)

pheatmap(mat =RNA_seq_matrix.log,border_color = F,cluster_rows = F)
pheatmap(mat =log2(RNA_seq_matrix + 1),border_color = F,cluster_rows = F)

##################################################################################################
# 根据Hela的数据，确定ribosome protein的表达量
##################################################################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")
RP_gene_table = read.table("~/menghw_HD/reference/ribosome_protein/rProtein_geneinfo.csv",header = T,sep = ";")
chrom_name_list = as.character(RP_gene_table$Chr_num)
chrom_name_list = sub("X","23",chrom_name_list)
chrom_name_list = sub("Y","24",chrom_name_list)
chrom_name_list = chrom.name.factor(as.integer(chrom_name_list))

hg19_gene_table = read.table("~/menghw_HD/reference/gene_gtf/hg19_annotation.gene.bed",header = T,sep = "\t")

RP_gene_table.fix = data.frame(chrom_name = chrom_name_list,
                               start = RP_gene_table$Start,
                               end = RP_gene_table$End,
                               gene_name = RP_gene_table$GeneName,
                               gene_type = RP_gene_table$Type)

hela_gene_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/Hela_gene_exp_all.table",header = T,sep = "\t")

RP_gene_exp_table = NULL
for(row.index in c(1:nrow(RP_gene_table.fix))){
  gene_name = as.character(RP_gene_table.fix$gene_name[row.index])
  table_index = which(hela_gene_exp_table$gene_name == gene_name)
  if(length(table_index) == 1){
    RP_gene_table.row = cbind(RP_gene_table.fix[row.index,],hela_gene_exp_table[row.index,c(7:10,12:19)])
    RP_gene_exp_table = rbind(RP_gene_exp_table,RP_gene_table.row)
  }
}

# 将结果保存成table
write.table(RP_gene_exp_table,file="~/menghw_HD/data_table/hela_ribosome_protein_gene_exp.table",col.names = T,row.names = F,sep = "\t",quote = F)

# 画RP gene expression的heatmap
RP_gene_exp_table = read.table("~/menghw_HD/data_table/hela_ribosome_protein_gene_exp.table",header = T,sep = "\t")
RP_gene_matrix = as.matrix(RP_gene_exp_table[,c(12:17)])
rownames(RP_gene_matrix) = as.character(RP_gene_exp_table$gene_name)

sd_vector = rep(0,nrow(RP_gene_exp_table))
for(row.index in c(1:nrow(RP_gene_matrix))){
  sd_vector[row.index] = sd(RP_gene_matrix[row.index,])
}
RP_gene_matrix.fix = RP_gene_matrix[sd_vector>0,]

pheatmap(RP_gene_matrix.fix,scale = "row",treeheight_row = F)


##################################################################################################
# 根据Hela的数据，确定House Keeping gene的表达量
##################################################################################################
rm(list=ls())
hela_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/Hela_gene_exp_all.table",header = T,sep = "\t")
hk_gene_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/HK_gene_list.bed",header = F,sep = "\t")
colnames(hk_gene_table) = c("chrom_name","start","end","gene_id","type","strand","gene_name")

hk_gene_table.fix = hela_exp_table[hela_exp_table$gene_id %in% hk_gene_table$gene_id,]
write.table(hk_gene_table.fix,file="~/menghw_HD/data_table/Hela-RNA-seq/HK_gene_exp.table",col.names = T,row.names = F,quote = F,sep = "\t")

library(pheatmap)
hk_gene_table.fix.matrix = as.matrix(hk_gene_table.fix[,c(14:19)])
colnames(hk_gene_table.fix.matrix) = colnames(hk_gene_table.fix[,c(14:19)])
rownames(hk_gene_table.fix.matrix) = hk_gene_table.fix$gene_name

pheatmap(log2(hk_gene_table.fix.matrix+1))
pheatmap(hk_gene_table.fix.matrix,scale = "row")


sd_vector = rep(0,nrow(hk_gene_table.fix.matrix))
for(row.index in c(1:nrow(hk_gene_table.fix.matrix))){
  sd_vector[row.index] = sd(hk_gene_table.fix.matrix[row.index,])
}
hk_gene_table.fix.matrix.fix = hk_gene_table.fix.matrix[sd_vector>0,]

pheatmap(hk_gene_table.fix.matrix.fix,scale = "row")

hk_gene_table.fix.significant = hk_gene_table.fix[hk_gene_table.fix$significant=="yes",]


RP_gene_table = hela_exp_table[unique(c(grep(pattern = "^RPS.*",hela_exp_table$gene_name) ,grep(pattern = "^RPL.*",hela_exp_table$gene_name))),]


##################################################################################################
# RP蛋白与NAT的关系
##################################################################################################
rm(list=ls())
load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")

hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
RP_gene_exp_table = read.table("~/menghw_HD/data_table/hela_ribosome_protein_gene_exp.table",header = T,sep = "\t")

chrom_name_list = chrom.name.factor(c(1:23))

overlap.table = overlap.gene(region_table = hela_NAT_table,RP_gene_exp_table)

gene_name_list = as.character(overlap.table$name.2)

RP_NAT_overlap = RP_gene_exp_table[RP_gene_exp_table$gene_name %in% gene_name_list,]


##################################################################################################
# house keeping gene in NAT
##################################################################################################
rm(list=ls())

load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")
hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
hk_gene_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/HK_gene_exp.table",header = T,sep = "\t")

NAT_hk_overlap = overlap(hela_NAT_table,hk_gene_table,level_vetor = chrom.name.factor(c(1:23)))

NAT_hk_gene_table = hk_gene_table[as.character(hk_gene_table$gene_name) %in% as.character(NAT_hk_overlap$name.2),]


##################################################################################################
# NAT gene table 
##################################################################################################
rm(list=ls())

load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")
hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
hela_gene_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/Hela_gene_exp_all.table",header = T,sep = "\t")

# function defination
overlap.gene <- function(region_table,gene_table,chrom_name_list=c(paste("chr",c(1:22),sep = ""),"chrX")){
  overlap.gene_table = NULL
  for(chrom_name in chrom_name_list){
    print(chrom_name)
    # 选择子集
    chrom_overlap.gene_table = NULL
    chrom_region_table = region_table[region_table[,1]==chrom_name,]
    chrom_gene_table = gene_table[gene_table[,1]==chrom_name,]
    for(chrom_region_index in c(1:nrow(chrom_region_table))){
      region_start = chrom_region_table[chrom_region_index,2]
      region_end = chrom_region_table[chrom_region_index,3]
      # length(which(!(chrom_gene_table[,2] > region_end | chrom_gene_table[,3] < region_start)))
      # region_filter.start = (chrom_gene_table[,2] >= region_start) & (chrom_gene_table[,2] <= region_end)
      # region_filter.end = (chrom_gene_table[,3] >= region_start) & (chrom_gene_table[,3] <= region_end)
      # region_filter = region_filter.start | region_filter.end
      # region_filter = region_filter.start
      region_filter = which(!(chrom_gene_table[,2] > region_end | chrom_gene_table[,3] < region_start))
      region_gene_table = chrom_gene_table[region_filter,]
      chrom_overlap.gene_table = rbind(chrom_overlap.gene_table,region_gene_table)
    }
    overlap.gene_table = rbind(overlap.gene_table,chrom_overlap.gene_table)  
  }
  
  return(overlap.gene_table)
}

NAT_gene_overlap_table = overlap.gene(region_table = hela_NAT_table,gene_table = hela_gene_exp_table)
NAT_gene_overlap_table.fix = NAT_gene_overlap_table[NAT_gene_overlap_table$gene_type=="protein_coding" & NAT_gene_overlap_table$test_state=="OK" ,]

NAT_gene_overlap_matrix = as.matrix(NAT_gene_overlap_table.fix[,c(14:19)])
colnames(NAT_gene_overlap_matrix) = colnames(NAT_gene_overlap_table.fix[,c(14:19)])
rownames(NAT_gene_overlap_matrix) = NAT_gene_overlap_table.fix$gene_name

sd_vector = rep(0,nrow(NAT_gene_overlap_matrix))
for(row.index in c(1:nrow(NAT_gene_overlap_matrix))){
  sd_vector[row.index] = sd(NAT_gene_overlap_matrix[row.index,])
}
NAT_gene_overlap_matrix.fix = NAT_gene_overlap_matrix[sd_vector>0,]

png(file="~/menghw_HD/R_image/Hela-hic-figure/NAIR_gene_heatmap.png",width = 500,height = 1600)
pheatmap(NAT_gene_overlap_matrix.fix,scale = "row",treeheight_row = F)
dev.off()

######################################################################################################
# NAIR gene changes v.s. other gene changes
######################################################################################################




######################################################################################################
# pheatmap and label
######################################################################################################

################################################
# HK gene with NAIR label
################################################


library(pheatmap)
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3 
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2 
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4 
colnames(test) = paste("Test", 1:10, sep = "") 
rownames(test) = paste("Gene", 1:20, sep = "")
# 设置每一列的注释
annotation_col = data.frame(CellType = factor(rep(c("CT1", "CT2"), 5)), Time = 1:5)
rownames(annotation_col) = paste("Test", 1:10, sep = "") # 设置每一行的注释
annotation_row = data.frame(GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6))))
rownames(annotation_row) = paste("Gene", 1:20, sep = "") # 设置注释的颜色
ann_colors = list(Time = c("white", "firebrick"),
                  CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
                  GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E"))
pheatmap(test,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors)






