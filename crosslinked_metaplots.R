#### crosslinked_metaplot.R
### 8/24/18
### Maria Mikedis mikedis@wi.mit.edu	
### run as: Rscript crosslinked_metaplots.R ${MOTIF}
### where MOTIF is a different motif sequence for each file (i.e., GUU, UUU, etc)
### each file has headers: nt.from.CL.site	no.motifs.CL	no.motifs.bkgd	no.CL.sites	no.bkgd.seq
### each files has rows ranging from -10 to 10


args <- commandArgs()
print(args)

motif <- args[6]  ### region is passed as argument #6; I see this with print(args)
library(ggplot2)

file_name = paste("CL_bkgd_", motif, "_counts" ,sep = "")
counts = read.table(file_name, header=TRUE, stringsAsFactors=TRUE)
counts = counts[!counts$nt.from.CL.site %in% c("n/a"),]
counts$nt.from.CL.site=as.numeric(levels(counts$nt.from.CL.site))[counts$nt.from.CL.site] ### keep order of X axis

counts$freq.CL = counts$no.motifs.CL/counts$no.CL.sites
counts$freq.bkgd = counts$no.motifs.bkgd/counts$no.bkgd.seq
counts$freq.motif = counts$freq.CL/counts$freq.bkgd


position = read.table(file_name)

output_file_name = paste("freq_", motif, ".pdf" ,sep = "")
title = paste(motif)

####################
## "Theme" for figures
####################

theme_mk <- theme(text=element_text(family="Helvetica",face="plain",size=6),
                  axis.text=element_text(size=6),
                  axis.title=element_text(size=6),
                  axis.line.y=element_blank(),
                  axis.line.x=element_line(size=0.25),
                  axis.ticks=element_line(size=0.25),
                  legend.title=element_blank(),
                  legend.background=element_rect(fill="transparent"),
                  legend.text=element_text(size=6), plot.title=element_text(size=6, face="plain"),
                  panel.background = element_blank(),
                  legend.key=element_blank())




ggplot(counts, aes(x = nt.from.CL.site, y = as.numeric(freq.motif), colour = "#80b1d3")) +
  geom_line(size=1) + 
  #geom_point(size = 4) + 
  theme_mk + 
  theme(#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0))) + 
  scale_color_manual(values=c("#80b1d3")) + ### I need this line and the original color designation in ggplot(aes()) to get color, for some reason
  ggtitle(title) +
  xlab("Distance from crosslinked nucleotide") +  
  ylab("Enrichment over background sequences") +  
  xlim(-10, 10)+
  ylim(0,20)
ggsave(output_file_name, scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)
