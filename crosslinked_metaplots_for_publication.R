#### crosslinked_metaplot_for_publication.R
### 8/24/18
### Maria Mikedis mikedis@wi.mit.edu	
### run as: Rscript crosslinked_metaplots_for_publication.R 
### plot for publication: one graph with motifs related to my major motif: GUU, UUU, UGUU(U/A)
## additional plot with other motifs: UGUU, GUU(U/A)
## one graph for motifs that were previously reported: GUUG (Zagore et al., 2018), GUUC (Maegawa et al., 2002 Genes to Cells; Reynolds et al., 20015 Human Mol Genetics), UUU[C/G]UUU) (Chen et al 2011 Genes and Dev Conti lab)


args <- commandArgs()
print(args)

library(ggplot2)

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




  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line.y=element_blank(), 
    axis.line.x = element_line(colour = "black"),
    text=element_text(size=20, family="Helvetica"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.title = element_text(family="Helvetica", size=20),
    legend.text = element_text(family="Helvetica", size=16))


## function to calculate motif frequency from file data
crosslinked_analysis <- function(motif) {
  file_name = paste("CL_bkgd_", motif, "_counts" ,sep = "")
  counts = read.table(file_name, header=TRUE, stringsAsFactors=TRUE)
  counts = counts[!counts$nt.from.CL.site %in% c("n/a"),]
  counts$nt.from.CL.site=as.numeric(levels(counts$nt.from.CL.site))[counts$nt.from.CL.site] ### keep order of X axis
  counts$freq.CL = counts$no.motifs.CL/counts$no.CL.sites
  counts$freq.bkgd = counts$no.motifs.bkgd/counts$no.bkgd.seq
  counts$freq.motif = counts$freq.CL/counts$freq.bkgd
  counts$motif = motif
  return(counts)
}

GUU_counts = crosslinked_analysis("GUU")
UUU_counts = crosslinked_analysis("UUU")
UGUU.UorA_counts = crosslinked_analysis("UGUU.UorA")
UGUU_counts = crosslinked_analysis("UGUU")
GUU.UorA_counts = crosslinked_analysis("GUU.UorA")
GUUG_counts = crosslinked_analysis("GUUG")
GUUC_counts = crosslinked_analysis("GUUC")
UUU.CorG.UUU_counts = crosslinked_analysis("UUU.CorG.UUU")


motifs_counts = rbind(UGUU.UorA_counts, UUU_counts, GUU_counts)
## order data
motifs_counts$motif = factor(motifs_counts$motif, c("UGUU.UorA", "GUU", "UUU"))


ggplot(motifs_counts, aes(x = nt.from.CL.site, y = freq.motif, colour = motif)) +
  geom_line(size=0.5) + 
  theme_mk + 
  scale_colour_manual(values=c("#386cb0", "#7fc97f","#fdc086"),  
            				  name="Motif",
                         breaks=c("UGUU.UorA", "GUU", "UUU"),
                         labels=c("UGUU(U/A)", "GUU", "UUU")) +  
  xlab("Distance from crosslinked nucleotide") +  
  ylab("Enrichment relative to background sequences") +  
  xlim(-10, 10)+
  ylim(0,25)
ggsave("metaplot_figure.pdf", scale = 1, width = 75, height=45, units = c("mm"), dpi = 300)

motifs_counts = rbind(UGUU.UorA_counts, UGUU_counts, GUU.UorA_counts)
## order data
motifs_counts$motif = factor(motifs_counts$motif, c("UGUU.UorA", "UGUU", "GUU.UorA"))

ggplot(motifs_counts, aes(x = nt.from.CL.site, y = freq.motif, colour = motif)) +
  geom_line(size=0.5) + 
  scale_colour_manual(values=c("#386cb0", "#7fc97f","#fdc086"),  
            				  name="Motif",
                         breaks=c("UGUU.UorA", "UGUU", "GUU.UorA"),
                         labels=c("UGUU(U/A)", "UGUU", "GUU(U/A)")) +  
  theme_mk + 
  xlab("Distance from crosslinked nucleotide") +  
  ylab("Enrichment relative to background sequences") +  
  xlim(-10, 10)+
  ylim(0,25)
ggsave("metaplot_figure2.pdf", scale = 1, width = 75, height=45, units = c("mm"), dpi = 300)


motifs_counts = rbind(UGUU.UorA_counts, GUUG_counts, GUUC_counts, UUU.CorG.UUU_counts)
## order data
motifs_counts$motif = factor(motifs_counts$motif, c("UGUU.UorA", "GUUG", "GUUC", "UUU.CorG.UUU"))

ggplot(motifs_counts, aes(x = nt.from.CL.site, y = freq.motif, colour = motif)) +
  geom_line(size=0.5) +  
  scale_colour_manual(values=c("#386cb0", "#7fc97f","#fdc086", "#beaed4"),  
            				  name="Motif",
                         breaks=c("UGUU.UorA", "GUUG", "GUUC",  "UUU.CorG.UUU"),
                         labels=c("UGUU(U/A)", "GUUG", "GUUC", "UUU(C/G)UUU")) +  
  theme_mk + 
  xlab("Distance from crosslinked nucleotide") +  
  ylab("Enrichment relative to background sequences") +  
  xlim(-10, 10)+
  ylim(0,25)
ggsave("metaplot_figure3.pdf", scale = 1, width = 75, height=45, units = c("mm"), dpi = 300)


