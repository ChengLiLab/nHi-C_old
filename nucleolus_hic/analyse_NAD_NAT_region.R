#############################################################################################################
# analysis of NAD/NAT with B compartment
#############################################################################################################
rm(list=ls())

##################################################################
# loading data 
##################################################################
# sequencing depth
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")

# A/B compartment with 100Kb binsize
hela_AB_table = read.table("~/menghw_HD/data_table/Hela-genome-hic_AB_all_100000_ice.table",header = T,sep = "\t")

# NAD/NAT region 
hela_NAD_table = read.table("~/menghw_HD/data_table/Hela_NAD.bed",header = F,sep = "\t")
hela_NAT_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")

# Hi-C column Sum 20Kb binsize
hic_colsum_table = read.table("~/menghw_HD/data_table/Hela_hic_colsum.table",header = T,sep = "\t")

# function 
load(file="~/menghw_HD/R_code/my_function/MyPlot_v06.RData")


##################################################################
# B compartment overlap(有多少NAIR在B中)
##################################################################
chrom_name_list = c(paste("chr",c(1:22),sep = ""),"chrX")

NAD_AB_overlap_table = overlap(hela_NAD_table,hela_AB_table,level_vetor = chrom_name_list)
NAD_AB_overlap_table.B  = NAD_AB_overlap_table[NAD_AB_overlap_table$name.2=="B",]
sum(NAD_AB_overlap_table.B$end - NAD_AB_overlap_table.B$start) / sum(hela_NAD_table$V3 - hela_NAD_table$V2)

NAD_B_overlap.ratio = rep(0,23)
for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_NAD_table = hela_NAD_table[hela_NAD_table$V1==chrom_name,]
  chrom_B_overlap_table = NAD_AB_overlap_table[NAD_AB_overlap_table$name.2=="B" & NAD_AB_overlap_table$chrom_name==chrom_name,]
  NAD_B_overlap.ratio[chrom_index] = sum(chrom_B_overlap_table$end -  chrom_B_overlap_table$start)/sum(chrom_NAD_table[,3] - chrom_NAD_table[,2])
}

barplot(NAD_B_overlap.ratio,names.arg = c(1:23),col="#7080D7")
length(which(hela_AB_table$info=="B")) / nrow(hela_AB_table)
abline(h=0.481,lwd=3,lty=2)

hela_NAD_table[hela_NAD_table$V1=="chr4",]
hela_NAD_table[hela_NAD_table$V1=="chr13",]


NAT_AB_overlap_table = overlap(hela_NAT_table,hela_AB_table,level_vetor = chrom_name_list)
NAT_AB_overlap_table.B  = NAT_AB_overlap_table[NAT_AB_overlap_table$name.2=="B",]
sum(NAT_AB_overlap_table.B$end - NAT_AB_overlap_table.B$start) / sum(hela_NAT_table$V3 - hela_NAT_table$V2)

NAT_B_overlap.ratio = rep(0,23)
for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_NAT_table = hela_NAT_table[hela_NAT_table$V1==chrom_name,]
  chrom_B_overlap_table = NAT_AB_overlap_table[NAT_AB_overlap_table$name.2=="B" & NAT_AB_overlap_table$chrom_name==chrom_name,]
  NAT_B_overlap.ratio[chrom_index] = sum(chrom_B_overlap_table$end -  chrom_B_overlap_table$start)/sum(chrom_NAT_table[,3] - chrom_NAT_table[,2])
}


######################################################
# 计算每个NAT中的B ratio  dotplot
######################################################
NAT_B_ratio = rep(0,23)
hela_NAT_table.fix = hela_NAT_table[(hela_NAT_table$V3-hela_NAT_table$V2) >= 1000e3,]

for(chrom_index in c(1:23)){
  print(chrom_index)
  chrom_name = chrom.name(chrom_index)
  chrom_NAT_table = hela_NAT_table.fix[hela_NAT_table.fix$V1==chrom_name,]
  chrom_NAT.overlap = overlap(chrom_NAT_table,hela_AB_table,level_vetor = chrom_name)
  chrom_NAT.overlap.B = chrom_NAT.overlap[chrom_NAT.overlap$name.2=="B",]
  NAT_B_ratio[chrom_index] = sum(chrom_NAT.overlap.B$end - chrom_NAT.overlap.B$start) / sum(chrom_NAT_table[,3] - chrom_NAT_table[,2])
}

NAT_B_ratio.table = data.frame(chrom_name = c(paste("chr",c(1:22),sep = ""),"chrX"),
                               NAT_B_ratio = NAT_B_ratio)

## 按照染色体从小到大排序
HG19_LEN=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
           135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,
           51304566,155270560,59373566,42999)


NAT_B_ratio.table$chrom_name = factor(x = NAT_B_ratio.table$chrom_name,levels = chrom.name.factor(order(HG19_LEN[1:23])))

require(ggplot2)
# ggplot(hela_NAT_table.fix,aes(x=chrom_name,y=NAT_B_ratio)) + geom_boxplot() + geom_jitter(position=position_jitter(0.2),col="#FF7304") + theme(legend.position ="none") + theme_gray()
ggplot(NAT_B_ratio.table,aes(x=chrom_name,y=NAT_B_ratio))  + geom_jitter(col="#FF7304")


######################################################
# 计算每条染色体NAT中的B的ratio boxplot
######################################################

NAT_B_ratio = rep(0,nrow(hela_NAT_table))
for(row.index in c(1:nrow(hela_NAT_table))){
  print(row.index)
  row.region = hela_NAT_table[row.index,]
  row.chrom_name = as.character(row.region[1,1])
  row.overlap = overlap(row.region,hela_AB_table,level_vetor = row.chrom_name)
  row.overlap.B = row.overlap[row.overlap$name.2=="B",]
  row.ratio.B = sum(row.overlap.B$end - row.overlap.B$start) / sum(row.region[,3] - row.region[,2])
  NAT_B_ratio[row.index] = row.ratio.B
}
hela_NAT_table.fix = cbind(hela_NAT_table,NAT_B_ratio)
# write.table(hela_NAT_table.fix,file = "~/menghw_HD/data_table/Hela_NAT_B_ratio.bed",col.names = F,row.names = F,quote = F,sep = '\t')
colnames(hela_NAT_table.fix) = c("chrom_name","start","end","peak_index","value","NAT_B_ratio")
# hela_NAT_table.fix$chrom_name = factor(x = hela_NAT_table.fix$chrom_name,levels = c(paste("chr",c(1:22),sep = ""),"chrX"))


hela_NAT_table.fix.raw = hela_NAT_table.fix
hela_NAT_table.fix = hela_NAT_table.fix.raw[(hela_NAT_table$V3-hela_NAT_table$V2) >= 1000e3,]

## 按照median从小到大排序
chrom_NAT_median = rep(1,23)
for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_NAT_median[chrom_index] = median(hela_NAT_table.fix$NAT_B_ratio[hela_NAT_table.fix$chrom_name==chrom_name])
}
hela_NAT_table.fix$chrom_name = factor(x = hela_NAT_table.fix$chrom_name,levels = chrom.name.factor(order(chrom_NAT_median)))


## 按照染色体从小到大排序
HG19_LEN=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
           135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,
           51304566,155270560,59373566,42999)


hela_NAT_table.fix$chrom_name = factor(x = hela_NAT_table.fix$chrom_name,levels = chrom.name.factor(order(HG19_LEN[1:23])))


require(ggplot2)
ggplot(hela_NAT_table.fix,aes(x=chrom_name,y=NAT_B_ratio)) + geom_boxplot() + geom_jitter(position=position_jitter(0.2),col="#FF7304") + theme(legend.position ="none") + theme_gray()

ggplot(hela_NAT_table.fix,aes(x=chrom_name,y=NAT_B_ratio)) + geom_boxplot()  + theme(legend.position ="none") + theme_gray() 




######################################################
# 计算每个NAD中的B ratio
######################################################
NAD_B_ratio = rep(0,nrow(hela_NAD_table))
for(row.index in c(1:nrow(hela_NAD_table))){
  print(row.index)
  row.region = hela_NAD_table[row.index,]
  row.chrom_name = as.character(row.region[1,1])
  row.overlap = overlap(row.region,hela_AB_table,level_vetor = row.chrom_name)
  row.overlap.B = row.overlap[row.overlap$name.2=="B",]
  row.ratio.B = sum(row.overlap.B$end - row.overlap.B$start) / sum(row.region[,3] - row.region[,2])
  NAD_B_ratio[row.index] = row.ratio.B
}
hela_NAD_table.fix = cbind(hela_NAD_table,NAD_B_ratio)
write.table(hela_NAD_table.fix,file = "~/menghw_HD/data_table/Hela_NAD_B_ratio.bed",col.names = F,row.names = F,quote = F,sep = '\t')
colnames(hela_NAD_table.fix) = c("chrom_name","start","end","peak_index","value","NAD_B_ratio")
# hela_NAD_table.fix$chrom_name = factor(x = hela_NAD_table.fix$chrom_name,levels = c(paste("chr",c(1:22),sep = ""),"chrX"))

## 按照median从小到大排序
chrom_NAD_median = rep(1,23)
for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_NAD_median[chrom_index] = median(hela_NAD_table.fix$NAD_B_ratio[hela_NAD_table.fix$chrom_name==chrom_name])
}
hela_NAD_table.fix$chrom_name = factor(x = hela_NAD_table.fix$chrom_name,levels = chrom.name.factor(order(chrom_NAD_median)))

require(ggplot2)
ggplot(hela_NAD_table.fix,aes(x=chrom_name,y=NAD_B_ratio)) + geom_boxplot() + geom_jitter(position=position_jitter(0.2),col="#5C5393") + theme(legend.position ="none") + theme_gray()


chrom.name.factor(c(2,23,35,12))

barplot(NAT_B_overlap.ratio,names.arg = c(1:23),col = "#FF7304",ylim=c(0,1))
length(which(hela_AB_table$info=="B")) / nrow(hela_AB_table)
abline(h=0.481,lwd=3,lty=2)


##################################################################
# B compartment overlap(NAIR能解释多少B compartment)
##################################################################
NAT_B_ratio = rep(0,23)

hela_NAT_table.fix = hela_NAT_table
# hela_NAT_table.fix = hela_NAT_table[(hela_NAT_table$V3-hela_NAT_table$V2) >= 1000e3,]

NAT_B_ratio = rep(0,23)

for(chrom_index in c(1:23)){
  print(chrom_index)
  chrom_name = chrom.name(chrom_index)
  chrom_NAT_table = hela_NAT_table.fix[hela_NAT_table.fix$V1==chrom_name,]
  chrom_NAT.overlap = overlap(chrom_NAT_table,hela_AB_table,level_vetor = chrom_name)
  chrom_NAT.overlap.B = chrom_NAT.overlap[chrom_NAT.overlap$name.2=="B",]
  
  chrom_B_table = hela_AB_table[hela_AB_table$chrom_name==chrom_name & hela_AB_table$info=="B",]
  
  NAT_B_ratio[chrom_index] = sum(chrom_NAT.overlap.B$end - chrom_NAT.overlap.B$start) / sum(chrom_B_table[,3] - chrom_B_table[,2])
}

NAT_B_ratio.table = data.frame(chrom_name = c(paste("chr",c(1:22),sep = ""),"chrX"),
                               NAT_B_ratio = NAT_B_ratio)

## 按照染色体从小到大排序
HG19_LEN=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
           135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,
           51304566,155270560,59373566,42999)


NAT_B_ratio.table$chrom_name = factor(x = NAT_B_ratio.table$chrom_name,levels = chrom.name.factor(order(HG19_LEN[1:23])))

require(ggplot2)
# ggplot(hela_NAT_table.fix,aes(x=chrom_name,y=NAT_B_ratio)) + geom_boxplot() + geom_jitter(position=position_jitter(0.2),col="#FF7304") + theme(legend.position ="none") + theme_gray()
ggplot(NAT_B_ratio.table,aes(x=chrom_name,y=NAT_B_ratio))  + geom_bar(fill="#FF7304",stat = "identity") 



##################################################################
# signal strength
##################################################################
hela_NAT_table.fix = hela_NAT_table
hela_NAT_table.fix = cbind(hela_NAT_table.fix,rep(0,nrow(hela_NAT_table.fix)))
colnames(hela_NAT_table.fix) = c("chrom_name","start","end","peak_index","value","ratio")

# for(row.index in c(1:nrow(hela_NAT_table))){
  
for(row.index in c(1:nrow(hela_NAT_table))){  
  print(row.index)
  table.row = hela_NAT_table[row.index,]
  chrom_name = as.character(table.row[1,1])
  region_start = as.integer(table.row[1,2])
  region_end = as.integer(table.row[1,3])
  
  select_vector = hic_colsum_table$chrom_name == as.character(chrom_name) & hic_colsum_table$region_start >= region_start & hic_colsum_table$region_end <= region_end
  n_hic_colsum_vector = hic_colsum_table$n_hic_colsum[select_vector]
  g_hic_colsum_vector = hic_colsum_table$g_hic_colsum[select_vector]
  
  ratio = (sum(n_hic_colsum_vector) / hic_count.hela_n_hic.total) / (sum(g_hic_colsum_vector) / hic_count.hela_g_hic.total)
  
  hela_NAT_table.fix$ratio[row.index] = ratio
}


hela_NAT_table.fix = cbind(hela_NAT_table.fix,rep(0,nrow(hela_NAT_table.fix)))
colnames(hela_NAT_table.fix) = c("chrom_name","start","end","peak_index","value","ratio","index")
for(row.index in c(1:nrow(hela_NAT_table.fix))){
  print(row.index)
  table.row = hela_NAT_table[row.index,]
  chrom_name = as.character(table.row[1,1])
  if(chrom_name == "rDNA"){
    chrom_index = 25
  }else if(chrom_name == "chrY"){
    chrom_index = 24
  }else if(chrom_name == "chrX"){
    chrom_index = 23
  }else{
    chrom_index = as.integer(sub("chr","",chrom_name))
  }
  
  hela_NAT_table.fix[row.index,7] = chrom_index

}

hela_NAT_table.fix = hela_NAT_table.fix[hela_NAT_table.fix$ratio<=4.5,]
hela_NAT_table.fix = hela_NAT_table.fix[(hela_NAT_table.fix$end - hela_NAT_table.fix$start) >= 1000e3,]

## 按照median从小到大排序
chrom_NAT_median = rep(1,23)
for(chrom_index in c(1:23)){
  chrom_name = chrom.name(chrom_index)
  chrom_NAT_median[chrom_index] = median(hela_NAT_table.fix$ratio[hela_NAT_table.fix$chrom_name==chrom_name])
}
hela_NAT_table.fix$chrom_name = factor(x = hela_NAT_table.fix$chrom_name,levels = chrom.name.factor(order(chrom_NAT_median)))



require(ggplot2)
ggplot(hela_NAT_table.fix,aes(x=chrom_name,y=ratio)) + geom_boxplot() + geom_jitter(position=position_jitter(0.2),col="#5C5393") + theme(legend.position ="none") + theme_gray()


boxplot(ratio ~ index, data = hela_NAT_table.fix, col = "lightgray",ylim=c(1.4,4))


sum(NAT_AB_overlap_table.B$end - NAT_AB_overlap_table.B$start) / nrow(hela_AB_table[hela_AB_table$info=="B",]) / 100e3

##################################################################
# NAT B-compatment ratio
##################################################################
hela_NAT_table.fix = hela_NAT_table
hela_NAT_table.fix = cbind(hela_NAT_table.fix,rep(0,nrow(hela_NAT_table.fix)))
colnames(hela_NAT_table.fix) = c("chrom_name","start","end","peak_index","value","ratio")

# for(row.index in c(1:nrow(hela_NAT_table))){

for(row.index in c(1:nrow(hela_NAT_table))){  
  print(row.index)
  table.row = hela_NAT_table[row.index,]
  chrom_name = as.character(table.row[1,1])
  region_start = as.integer(table.row[1,2])
  region_end = as.integer(table.row[1,3])
  
  select_vector = hic_colsum_table$chrom_name == as.character(chrom_name) & hic_colsum_table$region_start >= region_start & hic_colsum_table$region_end <= region_end
  n_hic_colsum_vector = hic_colsum_table$n_hic_colsum[select_vector]
  g_hic_colsum_vector = hic_colsum_table$g_hic_colsum[select_vector]
  
  ratio = (sum(n_hic_colsum_vector) / hic_count.hela_n_hic.total) / (sum(g_hic_colsum_vector) / hic_count.hela_g_hic.total)
  
  hela_NAT_table.fix$ratio[row.index] = ratio
}


hela_NAT_table.fix = cbind(hela_NAT_table.fix,rep(0,nrow(hela_NAT_table.fix)))
colnames(hela_NAT_table.fix) = c("chrom_name","start","end","peak_index","value","ratio","index")
for(row.index in c(1:nrow(hela_NAT_table.fix))){
  print(row.index)
  table.row = hela_NAT_table[row.index,]
  chrom_name = as.character(table.row[1,1])
  if(chrom_name == "rDNA"){
    chrom_index = 25
  }else if(chrom_name == "chrY"){
    chrom_index = 24
  }else if(chrom_name == "chrX"){
    chrom_index = 23
  }else{
    chrom_index = as.integer(sub("chr","",chrom_name))
  }
  
  hela_NAT_table.fix[row.index,7] = chrom_index
  
}


boxplot(ratio ~ index, data = hela_NAT_table.fix, col = "lightgray",ylim=c(1.4,4))












