##################################################################################################
# rmdup 去掉gene 名重复的行
##################################################################################################
rm(list=ls())

hela_gene_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/Hela_gene_exp_all.table",header = T,sep = "\t")
unique_state = rep(T,nrow(hela_gene_exp_table))
hela_gene_exp_table = cbind(hela_gene_exp_table,unique_state)

for(row.index in c(1:nrow(hela_gene_exp_table))){
  gene_id = hela_gene_exp_table$gene_id[row.index]
  if(length(which(hela_gene_exp_table$gene_id == as.character(gene_id))) > 1){
    hela_gene_exp_table$unique_state[row.index] = FALSE
  }
}

hela_gene_exp_table.unique = hela_gene_exp_table[hela_gene_exp_table$unique_state==T,]
write.table(hela_gene_exp_table.unique,file = "~/menghw_HD/data_table/Hela-RNA-seq/Hela_gene_exp_all_rmdup.table",col.names = T,row.names = F,quote = F,sep = "\t")

##################################################################################################
# 查找
##################################################################################################


##################################################################################################
# overlap function
##################################################################################################
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


##################################################################################################
# NAT gene table 
##################################################################################################
rm(list=ls())

# load(file="~/menghw_HD/R_code/my_function/MyPlot_v07.RData")
load(file="~/menghw_HD/R_code/my_function/overlap_gene.RData")
hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
hela_gene_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/Hela_gene_exp_all_rmdup.table",header = T,sep = "\t")

NAT_gene_overlap_table = overlap.gene(region_table = hela_NAT_table,gene_table = hela_gene_exp_table)
NAT_gene_overlap_table.fix = NAT_gene_overlap_table[NAT_gene_overlap_table$gene_type=="protein_coding" & NAT_gene_overlap_table$test_state=="OK" ,]
# write.table(NAT_gene_overlap_table.fix,file = "~/menghw_HD/data_table/Hela-RNA-seq/Hela_NAT_gene_exp.table",col.names = T,row.names = F,sep = "\t",quote = F)

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

##################################################################################################
# house keeping gene changes 
##################################################################################################
rm(list=ls())

hela_NAT_gene_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/Hela_NAT_gene_exp.table",header = T,sep = "\t")
hk_gene_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/HK_gene_exp.table",header = T,sep = "\t")

unique_state = rep(T,nrow(hk_gene_exp_table))
hk_gene_exp_table = cbind(hk_gene_exp_table,unique_state)

for(row.index in c(1:nrow(hk_gene_exp_table))){
  gene_id = hk_gene_exp_table$gene_id[row.index]
  if(length(which(hk_gene_exp_table$gene_id == as.character(gene_id))) > 1){
    hk_gene_exp_table$unique_state[row.index] = FALSE
  }
}

hk_gene_exp_table.unique = hk_gene_exp_table[hk_gene_exp_table$unique_state==T,]


hk_gene_exp_matrix = as.matrix(hk_gene_exp_table.unique[,c(14:19)])
colnames(hk_gene_exp_matrix) = colnames(hk_gene_exp_table.unique[,c(14:19)])
rownames(hk_gene_exp_matrix) = hk_gene_exp_table.unique$gene_name

sd_vector = rep(0,nrow(hk_gene_exp_matrix))
for(row.index in c(1:nrow(hk_gene_exp_matrix))){
  sd_vector[row.index] = sd(hk_gene_exp_matrix[row.index,])
}
hk_gene_exp_matrix.fix = hk_gene_exp_matrix[sd_vector>0 & hk_gene_exp_table.unique$significant=="yes",]


gene_position = rep("Other",nrow(hk_gene_exp_matrix.fix))
gene_position[which(rownames(hk_gene_exp_matrix.fix) %in% as.character((hela_NAT_gene_exp_table$gene_name)))] = "NAIR"

annotation_row = data.frame(gene_position)
rownames(annotation_row) = rownames(hk_gene_exp_matrix.fix)

ann_colors = list(gene_position = c(NAIR="#FF7304",Other="#5CCDC9"))
pheatmap(hk_gene_exp_matrix.fix,scale = "row",treeheight_row = F,annotation_row = annotation_row,annotation_colors = ann_colors)

##################################################################################################
# ribosome protein expression
##################################################################################################
rm(list=ls())
# load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")

hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")
RP_gene_exp_table = read.table("~/menghw_HD/data_table/hela_ribosome_protein_gene_exp.table",header = T,sep = "\t")

RP_NAT_overlap = overlap.gene(region_table = hela_NAT_table,RP_gene_exp_table)
gene_name_list = as.character(RP_NAT_overlap$gene_name)

RP_gene_matrix = as.matrix(RP_gene_exp_table[,12:17])
rownames(RP_gene_matrix) = RP_gene_exp_table$gene_name

sd_vector = rep(0,nrow(RP_gene_matrix))
for(row.index in c(1:nrow(RP_gene_matrix))){
  sd_vector[row.index] = sd(RP_gene_matrix[row.index,])
}
RP_gene_matrix.fix = RP_gene_matrix[sd_vector>0,]
RP_gene_exp_table.fix = RP_gene_exp_table[sd_vector>0,]

gene_position = rep("Other",nrow(RP_gene_matrix.fix))
gene_position[RP_gene_exp_table.fix$gene_name %in% gene_name_list] = "NAIR"

annote_row = data.frame(gene_position)
rownames(annote_row) = rownames(RP_gene_matrix.fix)

annote_col = list(gene_position=c(NAIR="#FF7304",Other="#5CCDC9"))
pheatmap(RP_gene_matrix.fix,scale = "row",treeheight_row = F,annotation_row = annote_row,annotation_colors = annote_col)
