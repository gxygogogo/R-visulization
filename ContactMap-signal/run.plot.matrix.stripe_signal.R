library(RColorBrewer)
library(ggplot2)
background<-(
  theme_bw() +
    theme(panel.background = element_blank()) +
    theme(panel.border = element_rect(size=1.5, colour="#939393", fill=NA)) +
    theme(panel.grid = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(axis.title = element_blank()) +
    theme(axis.text = element_blank()) +
    theme(legend.position="none")
)

plot_compared_contact_map <- function(tab_signal, tab_rawMatrix, reg, res, point_size = 0.6){
  #color_pal <- c("white",brewer.pal(name = "Reds",n=9)[1:6])
  color_pal = c("#F8FCF6", "#F3FBEF", "#EEF9E9", "#D6F0CB", "#C7E9C0FF", "#A1D99BFF", "#74C476FF", "#41AB5DFF", "#238B45FF", "#006D2CFF", "#00441BFF")
  
  tab_signal$bin1 <- tab_signal$bin1
  tab_signal$bin2 <- tab_signal$bin2
  
  tab_rawMatrix$bin1 <- tab_rawMatrix$bin1-1
  tab_rawMatrix$bin2 <- tab_rawMatrix$bin2-1
  
  tab_rawMatrix_down <- data.frame(chr=tab_rawMatrix$chr,
                                   bin1=tab_rawMatrix$bin2,
                                   bin2=tab_rawMatrix$bin1,
                                   delta=tab_rawMatrix$delta)
  
  tab_rawMatrix_up <- data.frame(chr=tab_rawMatrix$chr,
                                 bin1=tab_rawMatrix$bin1,
                                 bin2=tab_rawMatrix$bin2,
                                 delta=tab_rawMatrix$delta)
  
  reg=gsub(",", "", reg)
  chr=strsplit(reg, split="[_:-]")[[1]][1]
  start=as.numeric(strsplit(reg, split="[_:-]")[[1]][2])
  end=as.numeric(strsplit(reg, split="[_:-]")[[1]][3])
  start = start - start%%res
  end = end + res - end%%res
  x_min=start
  y_min=start
  x_max=end
  y_max=end
  reg_frame_df <- data.frame(chr=rep(chr, 4),
                             bin1=c(x_min, x_min, x_max, x_max),
                             bin2=c(y_min, y_max, y_min, y_max),
                             delta=c(0,0,0,0))
  
  dfplot <- rbind(reg_frame_df,tab_rawMatrix_up,tab_rawMatrix_down)
  
  
  dfplot$delta[dfplot$delta > 8] <- 8
  ggplot(data=dfplot) + geom_tile(aes(x=bin2, y=-bin1, fill=delta, colour=delta)) +
    geom_abline(slope=-1, intercept=0, linetype=2, size=0.4, colour="#bdbdbd") +
    scale_fill_gradientn(colours=color_pal, na.value = "white") +
    scale_color_gradientn(colours=color_pal, na.value = "white") +
    scale_x_continuous(expand = c(0.005,0.005)) +
    scale_y_continuous(expand = c(0.005,0.005)) +
    background +
    geom_point(tab_signal[which(tab_signal$delta=="Stripe"),],mapping = aes(bin2, -bin1),shape=22,colour="red",fill="red", size = point_size, stroke = 0) +
    geom_point(tab_signal[which(tab_signal$delta=="Reseting"),],mapping = aes(bin2, -bin1),shape=22,colour="grey",fill="grey", size = point_size, stroke = 0) +
    geom_point(tab_signal[which(tab_signal$delta=="CTCF_surrounding"),],mapping = aes(bin2, -bin1),shape=22,colour="blue",fill="blue", size = point_size, stroke = 0)
  
}


args<-commandArgs(T)
reg <- args[1]
point_size <- as.numeric(args[2])
# reg <- "chr1_165013116_166618115"
res <- 5000
setwd(paste("/data/tang/cohesin_project/ChIA_PET_FINAL/Final.for.Downstream.Analysis/run.call.signals/summary_of_signal/categorized.signals.for.visualization/", reg, sep=""))

tab_signal <- read.table(file=paste0("WT_SMC1A.signal.",reg,".bedpe"), col.names=c('chr', 'bin1', 'bin2', 'delta'))
tab_rawMatrix <- read.table(file=paste0("WT_SMC1A.rawMatrix.",reg,".bedpe"), col.names=c('chr', 'bin1', 'bin2', 'delta'))
WT <- plot_compared_contact_map(tab_signal, tab_rawMatrix, reg, res, point_size)

tab_signal <- read.table(file=paste0("SA1KO_SMC1A.signal.",reg,".bedpe"), col.names=c('chr', 'bin1', 'bin2', 'delta'))
tab_rawMatrix <- read.table(file=paste0("SA1KO_SMC1A.rawMatrix.",reg,".bedpe"), col.names=c('chr', 'bin1', 'bin2', 'delta'))
SA1KO <- plot_compared_contact_map(tab_signal, tab_rawMatrix, reg, res, point_size)

tab_signal <- read.table(file=paste0("SA2KO_SMC1A.signal.",reg,".bedpe"), col.names=c('chr', 'bin1', 'bin2', 'delta'))
tab_rawMatrix <- read.table(file=paste0("SA2KO_SMC1A.rawMatrix.",reg,".bedpe"), col.names=c('chr', 'bin1', 'bin2', 'delta'))
SA2KO <- plot_compared_contact_map(tab_signal, tab_rawMatrix, reg, res, point_size)

pdf(paste0("signal.",reg,".pdf"),width = 15,height = 5)
cowplot::plot_grid(WT,SA1KO,SA2KO,ncol = 3)
dev.off()

