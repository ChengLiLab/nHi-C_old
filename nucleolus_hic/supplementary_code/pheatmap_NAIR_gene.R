#######################################################################################
# hela NAIR gene expression heatmap
#######################################################################################

############################################################
# function
############################################################
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

chrom.name.factor <- function(chrom_index_list){
  chrom_index_list = as.vector(chrom_index_list)
  chrom_name_list = rep("chrN",length(chrom_index_list))
  for(index in c(1:length(chrom_index_list))){
    chrom_index = chrom_index_list[index]
    if(chrom_index < 23){
      chrom_name = paste0("chr",chrom_index)
    }else if(chrom_index == 23){
      chrom_name = "chrX"
    }else if(chrom_index == 24){
      chrom_name = "chrY"
    }else if(chrom_index == 25){
      chrom_name = "rDNA"
    }else{
      chrom_name = NA
    }
    chrom_name_list[index] = chrom_name
  }
  return(chrom_name_list)
}

############################################################
# load data
############################################################
load(file="~/menghw_HD/R_code/my_function/overlap_gene.RData")
hela_gene_exp_table = read.table("~/menghw_HD/data_table/Hela-RNA-seq/Hela_gene_exp_all_rmdup.table",header = T,sep = "\t")
hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")

length(hela_gene_exp_table$gene_name)


chrom_index_list = c(1:23)
chrom_name_list = chrom.name.factor(chrom_index_list)

NAIR_gene_table = overlap.gene(hela_NAIR_table,hela_gene_exp_table,chrom_name_list = chrom_name_list)
NAIR_gene_table.filter = NAIR_gene_table[NAIR_gene_table$test_state=="OK" & NAIR_gene_table$gene_type=="protein_coding",]

require(pheatmap)
NAIR_gene_table.filter.matrix = as.matrix(NAIR_gene_table.filter[,c(14:19)])
rownames(NAIR_gene_table.filter.matrix) = NAIR_gene_table.filter$gene_name

pheatmap(log2(NAIR_gene_table.filter.matrix),scale = "row",cluster_rows = F,cluster_cols = T)
 






