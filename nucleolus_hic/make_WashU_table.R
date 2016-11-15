## function
rm(list=ls())
chrom.length <- function(chrom_index,GenomeVersion="hg19"){
  # 根据chr_index返回chr_len
  HG19_LEN=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
             135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,
             51304566,155270560,59373566,42999)
  
  HG38_LEN=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,
             135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,
             50818468,156040895,57227415,42999)
  if(GenomeVersion=="hg19"){
    return(HG19_LEN[chrom_index])
  }else if(GenomeVersion=="hg38"){
    return(HG38_LEN[chrom_index])
  }
}

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


## read table
hela_all_table = read.table("~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep = "\t")

hela_all_table.fix = hela_all_table
hela_all_table.fix$region_max = hela_all_table.fix$region_max - 1

# hela-g-seq
hela_g_seq_total_count = sum(hela_all_table.fix$g_count) / 1e6
binsize_len = 100e3 / 1e3
loci_count_vector = round((hela_all_table.fix$g_count / hela_g_seq_total_count / binsize_len * 2),4)
name_vector = rep("Hela-g-seq",nrow(hela_all_table.fix))
hela_g_seq_table.bed = data.frame(chrom_name = hela_all_table.fix$chr_index,
                                  chrom_start = hela_all_table.fix$region_min,
                                  chrom_end = hela_all_table.fix$region_max,
                                  value = loci_count_vector)

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_len = chrom.length(chrom_index)
  index = which(hela_g_seq_table.bed$chrom_name==chrom_name)[length(which(hela_g_seq_table.bed$chrom_name==chrom_name))]  
  hela_g_seq_table.bed$chrom_end[index] = chrom_len - 1
}

write.table(hela_g_seq_table.bed,file = "~/menghw_HD/data_table/WashU_File/Hela-g-seq_RPKM.bed",col.names = F,row.names = F,sep = "\t",quote = F)

# hela-n-seq
hela_n_seq_total_count = sum(hela_all_table.fix$n_count) / 1e6
binsize_len = 100e3 / 1e3
loci_count_vector = round((hela_all_table.fix$n_count / hela_n_seq_total_count / binsize_len * 2),4)

hela_n_seq_table.bed = data.frame(chrom_name = hela_all_table.fix$chr_index,
                                  chrom_start = hela_all_table.fix$region_min,
                                  chrom_end = hela_all_table.fix$region_max,
                                  value = loci_count_vector)

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_len = chrom.length(chrom_index)
  index = which(hela_n_seq_table.bed$chrom_name==chrom_name)[length(which(hela_n_seq_table.bed$chrom_name==chrom_name))]  
  hela_n_seq_table.bed$chrom_end[index] = chrom_len - 1
}

write.table(hela_n_seq_table.bed,file = "~/menghw_HD/data_table/WashU_File/Hela-n-seq_RPKM.bed",col.names = F,row.names = F,sep = "\t",quote = F)

chrom.length(16)



# hela-n-seq-ActD
hela_n_seq_ActD_total_count = sum(hela_all_table.fix$n_actd_count) / 1e6
binsize_len = 100e3 / 1e3
loci_count_vector = round((hela_all_table.fix$n_actd_count / hela_n_seq_ActD_total_count / binsize_len * 2),4)

hela_n_seq_ActD_table.bed = data.frame(chrom_name = hela_all_table.fix$chr_index,
                                  chrom_start = hela_all_table.fix$region_min,
                                  chrom_end = hela_all_table.fix$region_max,
                                  value = loci_count_vector)

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_len = chrom.length(chrom_index)
  index = which(hela_n_seq_ActD_table.bed$chrom_name==chrom_name)[length(which(hela_n_seq_ActD_table.bed$chrom_name==chrom_name))]  
  hela_n_seq_ActD_table.bed$chrom_end[index] = chrom_len - 1
}

write.table(hela_n_seq_ActD_table.bed,file = "~/menghw_HD/data_table/WashU_File/Hela-n-seq-ActD_RPKM.bed",col.names = F,row.names = F,sep = "\t",quote = F)



# hela-g-hic
hela_g_hic_total_count = sum(hela_all_table.fix$g_hic) / 1e6
binsize_len = 100e3 / 1e3
loci_count_vector = round((hela_all_table.fix$g_hic / hela_g_hic_total_count / binsize_len * 2),4)

hela_g_hic_table.bed = data.frame(chrom_name = hela_all_table.fix$chr_index,
                                  chrom_start = hela_all_table.fix$region_min,
                                  chrom_end = hela_all_table.fix$region_max,
                                  value = loci_count_vector)

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_len = chrom.length(chrom_index)
  index = which(hela_g_hic_table.bed$chrom_name==chrom_name)[length(which(hela_g_hic_table.bed$chrom_name==chrom_name))]  
  hela_g_hic_table.bed$chrom_end[index] = chrom_len - 1
}

write.table(hela_g_hic_table.bed,file = "~/menghw_HD/data_table/WashU_File/Hela-g-hic_RPKM.bed",col.names = F,row.names = F,sep = "\t",quote = F)



# hela-n-hic
hela_n_hic_total_count = sum(hela_all_table.fix$n_hic) / 1e6
binsize_len = 100e3 / 1e3
loci_count_vector = round((hela_all_table.fix$n_hic / hela_n_hic_total_count / binsize_len * 2),4)

hela_n_hic_table.bed = data.frame(chrom_name = hela_all_table.fix$chr_index,
                                  chrom_start = hela_all_table.fix$region_min,
                                  chrom_end = hela_all_table.fix$region_max,
                                  value = loci_count_vector)

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_len = chrom.length(chrom_index)
  index = which(hela_n_hic_table.bed$chrom_name==chrom_name)[length(which(hela_n_hic_table.bed$chrom_name==chrom_name))]  
  hela_n_hic_table.bed$chrom_end[index] = chrom_len - 1
}

write.table(hela_n_hic_table.bed,file = "~/menghw_HD/data_table/WashU_File/Hela-n-hic_RPKM.bed",col.names = F,row.names = F,sep = "\t",quote = F)



# hela-n-hic-ActD
hela_n_hic_ActD_total_count = sum(hela_all_table.fix$n_hic_ActD) / 1e6
binsize_len = 100e3 / 1e3
loci_count_vector = round((hela_all_table.fix$n_hic_ActD / hela_n_hic_ActD_total_count/ binsize_len * 2),4)

hela_n_hic_ActD_table.bed = data.frame(chrom_name = hela_all_table.fix$chr_index,
                                  chrom_start = hela_all_table.fix$region_min,
                                  chrom_end = hela_all_table.fix$region_max,
                                  value = loci_count_vector)

for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_len = chrom.length(chrom_index)
  index = which(hela_n_hic_ActD_table.bed$chrom_name==chrom_name)[length(which(hela_n_hic_ActD_table.bed$chrom_name==chrom_name))]  
  hela_n_hic_ActD_table.bed$chrom_end[index] = chrom_len - 1
}

write.table(hela_n_hic_ActD_table.bed,file = "~/menghw_HD/data_table/WashU_File/Hela-n-hic-ActD_RPKM.bed",col.names = F,row.names = F,sep = "\t",quote = F)








