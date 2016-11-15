##################################################################################################
# compare two methods to call NAIR
##################################################################################################
rm(list=ls())

######################################################################
# No.1 use nHi-C / Hi-C column sum 
######################################################################
hela_NAIR_table = read.table("~/menghw_HD/data_table/Hela_NAT.bed",header = F,sep = "\t")

######################################################################
# No.2 use nHi-C / Seq column sum 
######################################################################
load(file="~/menghw_HD/data_table/R_data/sequence_depth.RData")
hela_all_table = read.table("~/menghw_HD/data_table/Hela_all_table_100000.txt",header = T,sep = "\t")

hela_g_seq.vector.fix = hela_all_table$g_count / sum(hela_all_table$g_count)
hela_g_hic.vector.fix = hela_all_table$g_hic / sum(hela_all_table$g_hic)
hela_n_hic.vector.fix = hela_all_table$n_hic / sum(hela_all_table$n_hic)

hela_NAIR.vector.1 = hela_n_hic.vector.fix / hela_g_hic.vector.fix
hela_NAIR.vector.1[hela_g_hic.vector.fix==0] = 0

hela_NAIR.vector.2 = hela_n_hic.vector.fix / hela_g_seq.vector.fix
hela_NAIR.vector.2[hela_g_seq.vector.fix==0] = 0

length(which(hela_NAIR.vector.2 > 1.4))

hela_NAIR_table.2 = hela_all_table[hela_NAIR.vector.2 > 1.7,]
write.table(hela_NAIR_table.2,file = "~/menghw_HD/data_table/raw_table/Hela_NAIR_2.bed",col.names = F,row.names = F,quote = F,sep = "\t")


######################################################################
# Plot ratio-ratio dot plot
######################################################################
par(mar=c(5,5,1,1))
plot(hela_NAIR.vector.1,hela_NAIR.vector.2,col=rgb(0,0,1,0.1),frame.plot=F,xaxt="n",yaxt="n",xlab="",ylab="",pch=19,xlim=c(0,10),ylim=c(0,10))
axis(side = 1,at=seq(0,10,2),labels =F,lwd=3,tck=-0.05)
axis(side = 2,at=seq(0,10,2),labels =F,lwd=3,tck=-0.05)



######################################################################
# calculate overlap
######################################################################
hela_NAIR_table.2.merge = read.table("~/menghw_HD/data_table/raw_table/Hela_NAIR_2_w1e5.bed",header = F,sep = "\t")
sum(hela_NAIR_table.2.merge$V3 - hela_NAIR_table.2.merge$V2)

load(file="~/menghw_HD/R_code/my_function/MyPlot_v08.RData")
hela_NAIR_overlap = overlap(hela_NAIR_table,hela_NAIR_table.2,level_vetor = c(paste("chr",c(1:22),sep = ""),"chrX"))
overlap_ratio = sum(hela_NAIR_overlap$end - hela_NAIR_overlap$start) / sum(hela_NAIR_table.2.merge$V3 - hela_NAIR_table.2.merge$V2) 







