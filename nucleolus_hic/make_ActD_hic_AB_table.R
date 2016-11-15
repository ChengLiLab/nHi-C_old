rm(list=ls())

require(HiTC)
require(GenomicRanges)
require(Matrix)

# 载入MyPlot 函数
load("~/menghw_HD/R_code/my_function/MyPlot_v07.RData")
image_path = "~/menghw_HD/R_image/Hela-hic-Rimage/Hela-hic_all_matrix_chr_%d_%d.RData"
png_path = "/home/menghw/menghw_HD/R_image/Hela-hic-figure/20160912/Hela-g-hic-ActD_confirm_AB_chr_%d_%d.png"
hela_ab_value_table = read.table("~/menghw_HD/data_table/Hela_AB_compart_value_fix.table",header = T,sep = "\t")

binsize = 100e3


# ## 开始循环
# for(chrom_index in c(7:12,15:22)){
#   load(sprintf(image_path,chrom_index,binsize))
#   chrom_len = chrom.length(chrom_index)
#   chrom_start = 0
#   chrom_end = chrom_len
#   chrom_name = chrom.name(chrom_index)
#   
# #   pca.method = ""
# #   if(is.null(pc.g_hic_ActD.norm)){
# #     pca.method = "mean"
# #     pc.g_hic.norm.fix <- pca.hic(g_hic_matrix.norm.obj, normPerExpected=TRUE, npc=1,method="loess")
# #     pc.g_hic_ActD.norm.fix <- pca.hic(g_hic_matrix.norm.obj, normPerExpected=TRUE, npc=1,method="loess")
# #     pc.g_hic_ActD.raw.fix <- pca.hic(g_hic_ActD_HiTCexp.obj, normPerExpected=TRUE, npc=1)
# #   }else{
# #     pca.method = "loess"
# #     pc.g_hic.norm.fix = pc.g_hic.norm
# #     pc.g_hic_ActD.norm.fix = pc.g_hic_ActD.norm
# #   }
#   
#   if(! is.null(pc.g_hic_ActD.norm)){
#   
#     png(sprintf(png_path,chrom_index,binsize),width = 1800,height = 2000)
#     
#     fig.facet <- layout(matrix(c(1:6),nrow = 3,byrow = FALSE),width = rep(10,6),heights = c(10,2,10,10,2,10))
#     layout.show(fig.facet)
#     
#     ## hela-g-hic
#     par(mar=c(1,1,1,1))
#     plot.matrix(g_hic_matrix.norm,bound.max=0.95)
#     
#     ## hela-g-hic A/B
#     par(mar=c(0,1,0,1),family="Arial",font=1)
#     ab_x = start(pc.g_hic.norm$PC1)
#     ab_y = score(pc.g_hic.norm$PC1)*(hela_ab_value_table$ice_coef[hela_ab_value_table$chrom_index==chrom_name])
#     plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
#     
#     ## hela-n-hic
#     par(mar=c(1,1,1,1),family="Arial",font=1)
#     plot.matrix(n_hic_matrix,bound.max = 0.95)
#     
#     ## hela-g-hic-ActD
#     par(mar=c(1,1,1,1),family="Arial",font=1)
#     plot.matrix(g_hic_ActD_matrix,bound.max = 0.95)
#     
#     ## hela-g-hic-ActD A/B
#     par(mar=c(0,1,0,1),family="Arial",font=1)
#     ab_x = start(pc.g_hic_ActD.norm$PC1)
#     ab_y = score(pc.g_hic_ActD.norm$PC1)
#     plot(ab_x,ab_y,type="h",xlab="", ylab="PC1vec", frame=F,xaxt="n",xlim=c(chrom_start,chrom_end),lwd=1,col=ifelse(ab_y>0,"orange","blue"))
#     
#     ## hela-n-hic
#     par(mar=c(1,1,1,1),family="Arial",font=1)
#     plot.matrix(n_hic_ActD_matrix,bound.max = 0.95)
#     
#     dev.off()
#   }
#   print(sprintf("The chr%d figure is done!",chrom_index))
# }


## 开始循环
# for(chrom_index in c(7:12,15:22)){
#   load(sprintf(image_path,chrom_index,binsize))
#   chrom_len = chrom.length(chrom_index)
#   chrom_start = 0
#   chrom_end = chrom_len
#   chrom_name = chrom.name(chrom_index)
  
  #   pca.method = ""
  #   if(is.null(pc.g_hic_ActD.norm)){
  #     pca.method = "mean"
  #     pc.g_hic.norm.fix <- pca.hic(g_hic_matrix.norm.obj, normPerExpected=TRUE, npc=1,method="loess")
  #     pc.g_hic_ActD.norm.fix <- pca.hic(g_hic_matrix.norm.obj, normPerExpected=TRUE, npc=1,method="loess")
  #     pc.g_hic_ActD.raw.fix <- pca.hic(g_hic_ActD_HiTCexp.obj, normPerExpected=TRUE, npc=1)
  #   }else{
  #     pca.method = "loess"
  #     pc.g_hic.norm.fix = pc.g_hic.norm
  #     pc.g_hic_ActD.norm.fix = pc.g_hic_ActD.norm
  #   }

hela_AB_ActD_table = NULL
for(chrom_index in c(7:12,15:22)){
  load(sprintf(image_path,chrom_index,binsize))
  chrom_name = chrom.name(chrom_index)
  
  # if(! is.null(pc.g_hic_ActD.norm)){
    
    region_start = start(pc.g_hic_ActD.norm$PC1)
    region_end = region_start + binsize
    region_value = score(pc.g_hic_ActD.norm$PC1) * hela_ab_value_table$ActD_ice_coef[hela_ab_value_table==chrom_name]
    region_type = rep("A",length(region_start))
    region_type[region_value < 0 ] = "B"
    
    chrom_AB_table = data.frame(chrom_name = rep(chrom_name,length(region_start)),
                                start = region_start,
                                end = region_end,
                                type = region_type,
                                PC1_value = region_value)
    
    hela_AB_ActD_table = rbind(hela_AB_ActD_table,chrom_AB_table)
  # }
  print(sprintf("The chr%d figure is done!",chrom_index))
}

write.table(hela_AB_ActD_table,file = "~/menghw_HD/data_table/Hela-g-hic-ActD_AB.table",col.names = T,sep = "\t",row.names = F,quote = F)











