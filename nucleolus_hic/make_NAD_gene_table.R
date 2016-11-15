rm(list=ls())
load("~/menghw_HD/R_code/my_function/MyPlot_v02.RData")

overlap.gene <- function(table_1,table_2,level_vetor=c("chr1")){
  # overlap 是针对排序过的table_1,table_2进行取overlap
  # 返回是1个overlap的data.frame
  overlap_table = data.frame(row.names = c("chrom_name","start","end","name","start.1","end.1","name.1","start.2","end.2","name.2","gene_type"))
  
  for(table.level in level_vetor){
    print(table.level)
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
      region_2.type = bed_table_2[index_2,5]
      
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
                                 name.2 = region_2.name,
                                 gene_type = region_2.type
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


NAD_table.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix.table",header = F,sep = "\t")
NAD_table.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-seq_narrow_fix_w1e4.table",header = T,sep = "\t")

# NAIR_table.fix = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix.table",header = F,sep = "\t")
# NAIR_table.merge = read.table("~/menghw_HD/data_table/fix_table/Hela-n-hic_broad_fix_w1e5.table",header = F,sep = "\t")

genome_gene_table = read.table("~/menghw_HD/reference/gene_gtf/hg19_annotation.gene.bed",header = T,sep = "\t")

chrom_name_vector = c(paste0("chr",c(1:22)),c("chrX"))

NAD_gene_table = overlap.gene(table_1 = NAD_table.merge,table_2 = genome_gene_table,level_vetor = chrom_name_vector )

# NAIR_gene_table = overlap.gene(table_1 = NAIR_table.merge,table_2 = genome_gene_table,level_vetor = chrom_name_vector )

write.table(NAD_gene_table,file = "~/menghw_HD/data_table/fix_table/NAD_gene_table.bed",col.names = T,row.names = F,sep = "\t",quote = F)

# write.table(NAIR_gene_table,file = "~/menghw_HD/data_table/fix_table/NAIR_gene_table.bed",col.names = T,row.names = F,sep = "\t",quote = F)






