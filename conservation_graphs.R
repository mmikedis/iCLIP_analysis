#### conservation_graphs.R
### 4/9/19
### Maria Mikedis mikedis@wi.mit.edu  
### run as: Rscript conservation_graphs.R REGION

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



##############
### phastCons analysis, crosslinked sites vs other 3' UTR nucleotides
##############



CL.phast = read.table("utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_phastCons.bed", header = FALSE, sep = "\t", stringsAsFactors =T)
bkgnd.phast = read.table("mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene_DAZL.bound_no.Dazl.bound.nts_phastCons.bed", header = FALSE, sep = "\t", stringsAsFactors =T)

# data is not formated as 1 nt per line; reformatting so that each line represents the score of 1 nt
CL.phast.reform = c()
i = 1
while (i <= nrow(CL.phast)) {
  m = rep(CL.phast$V4[i],each=(CL.phast$V3[i] - CL.phast$V2[i]))
  CL.phast.reform = c(CL.phast.reform,m)
  i=i+1
}

CL.phast.reform = as.data.frame(CL.phast.reform)
colnames(CL.phast.reform) = c("score")



bkgnd.phast.reform = c()
i = 1
while (i <= nrow(bkgnd.phast)) {
  m = rep(bkgnd.phast$V4[i],each=(bkgnd.phast$V3[i] - bkgnd.phast$V2[i]))
  bkgnd.phast.reform = c(bkgnd.phast.reform,m)
  i=i+1
}

bkgnd.phast.reform = as.data.frame(bkgnd.phast.reform)
colnames(bkgnd.phast.reform) = c("score")

CL.phast.reform$Dazl_target = paste("DAZL-bound nts")
bkgnd.phast.reform$Dazl_target = paste("unbound nts")


data = rbind(CL.phast.reform, bkgnd.phast.reform)

  ## order data
data$Dazl_target = factor(data$Dazl_target,c("DAZL-bound nts", "unbound nts"))



ggplot(data, aes(Dazl_target, score,fill=Dazl_target)) +
  geom_boxplot(outlier.size=0, outlier.stroke=0, lwd=0.2) +
  ylab("phastCons conservation score") +
  scale_fill_manual(values=c("#80b1d3", "gray")) +
  coord_cartesian(ylim=c(0, 1)) +
  scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75, 1)) + 
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))

ggsave("FIGURE_phastCons_CL.v.unbound_190409.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)


ggsave("FIGURE_phastCons_CL.v.unbound.v2_190409.pdf", scale = 2, width = 45, height=45, units = c("mm"), dpi = 300)

# Function to produce summary statistics (median and interquartile range)
data_summary <- function(x) {
   m <- median(x)
   ymin <- unname(quantile(x, 0.25))
   ymax <- unname(quantile(x, 0.75))
   return(c(y=m,ymin=ymin,ymax=ymax))
}

### violin plot
ggplot(data, aes(Dazl_target, score,fill=Dazl_target, colour=Dazl_target)) +
  geom_violin() +
  stat_summary(fun.data=data_summary, colour="black", size=0.005) +
  ylab("phastCons conservation score") +
  scale_fill_manual(values=c("#80b1d3", "gray")) +
   scale_colour_manual(values=c("#80b1d3", "gray")) +
  scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75, 1)) + 
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))



ggsave("FIGURE_phastCons_CL.v.unbound_violinplot_190409.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)


ggsave("FIGURE_phastCons_CL.v.unbound_violinplot.v2_190409.pdf", scale = 2, width = 45, height=45, units = c("mm"), dpi = 300)






## test for statistical difference between DAZL crosslinked and uncrosslinked nucleotides
stats= paste("Statistical difference of phastCons score between DAZL crosslinked nucleotides and uncrosslinked nucleotides")
stats = append(stats, capture.output(ks.test(CL.phast.reform$score, bkgnd.phast.reform$score)))
stats = append(stats, paste("\n\n"))
stats = append(stats, capture.output(wilcox.test(CL.phast.reform$score, bkgnd.phast.reform$score)))

writeLines(stats, "phastCons_CL.v.unbound_ks.test_wilcox.test_190408")

## get mean and median info for datasets
stats2= paste("Summary of CL nucleotide dataset")
stats2 = append(stats2, paste("\n"))
stats2 = append(stats2, capture.output(summary(CL.phast.reform)))
stats2 = append(stats2, paste("\n\n\n\n\n\n"))
stats2= append(stats2, paste("Summary of noncrosslinked nucleotide dataset"))
stats2 = append(stats2, paste("\n"))
stats2 = append(stats2, capture.output(summary(bkgnd.phast.reform)))

writeLines(stats2, "phastCons_CL.unbound_summary.stats_190408")




##############
### phyloP analysis, crosslinked sites vs other 3' UTR nucleotides
##############
CL.phylo = read.table("utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_phyloP.bed", header = FALSE, sep = "\t", stringsAsFactors =T)
bkgnd.phylo = read.table("mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene_DAZL.bound_no.Dazl.bound.nts_phyloP.bed", header = FALSE, sep = "\t", stringsAsFactors =T)




CL.phylo.reform = c()
i = 1
while (i <= nrow(CL.phylo)) {
  m = rep(CL.phylo$V4[i],each=(CL.phylo$V3[i] - CL.phylo$V2[i]))
  CL.phylo.reform = c(CL.phylo.reform,m)
  i=i+1
}
CL.phylo.reform = as.data.frame(CL.phylo.reform)
colnames(CL.phylo.reform) = c("score")



bkgnd.phylo.reform = c()
i = 1
while (i <= 2000) {
  m = rep(bkgnd.phylo$V4[i],each=(bkgnd.phylo$V3[i] - bkgnd.phylo$V2[i]))
  bkgnd.phylo.reform = c(bkgnd.phylo.reform,m)
  i=i+1
}

bkgnd.phylo.reform = as.data.frame(bkgnd.phylo.reform)
colnames(bkgnd.phylo.reform) = c("score")


CL.phylo.reform$Dazl_target = paste("DAZL-crosslinked nucleotides")
bkgnd.phylo.reform$Dazl_target = paste("uncrosslinked nucleotides")


data2 = rbind(CL.phylo.reform, bkgnd.phylo.reform)

  ## order data
data2$Dazl_target = factor(data2$Dazl_target,c("DAZL-crosslinked nucleotides", "uncrosslinked nucleotides"))



ggplot(data2, aes(Dazl_target, score,fill=Dazl_target)) +
  geom_boxplot(outlier.size=0, outlier.stroke=0, lwd=0.2) +
  ylab("phyloP conservation score") +
  scale_fill_manual(values=c("#80b1d3", "gray")) +
  coord_cartesian(ylim=c(-3.5, 6.5)) +
  scale_y_continuous(breaks=c(-3, -2, -1, 0,1,2,3,4,5,6)) + 
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))

ggsave("FIGURE_phyloP_CL.v.unbound_190409.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)


ggsave("FIGURE_phyloP_CL.v.unbound.v2_190409.pdf", scale = 2, width = 45, height=45, units = c("mm"), dpi = 300)




# Function to produce summary statistics (median and interquartile range)
data_summary <- function(x) {
   m <- median(x)
   ymin <- unname(quantile(x, 0.25))
   ymax <- unname(quantile(x, 0.75))
   return(c(y=m,ymin=ymin,ymax=ymax))
}

### violin plot
ggplot(data2, aes(Dazl_target, score,fill=Dazl_target, colour=Dazl_target)) +
  geom_violin() +
  stat_summary(fun.data=data_summary, colour="black", size=0.005) +
 ylab("phyloP conservation score") +
  scale_fill_manual(values=c("#80b1d3", "gray")) +
   scale_colour_manual(values=c("#80b1d3", "gray")) +
  coord_cartesian(ylim=c(-5, 7.5)) +
  scale_y_continuous(breaks=c(-5, -4, -3, -2, -1, 0,1,2,3,4,5,6, 7)) + 
  theme_mk + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"mm"),
        axis.title=element_text(margin=c(0,0,0,0)))



ggsave("FIGURE_phyloP_CL.v.unbound_violinplot_190409.pdf", scale = 1, width = 45, height=45, units = c("mm"), dpi = 300)


ggsave("FIGURE_phyloP_CL.v.unbound_violinplot.v2_190409.pdf", scale = 2, width = 45, height=45, units = c("mm"), dpi = 300)



## test for statistical difference between DAZL crosslinked and uncrosslinked nucleotides
stats= paste("Statistical difference of phyloP score between DAZL crosslinked nucleotides and uncrosslinked nucleotides")
stats = append(stats, capture.output(ks.test(CL.phylo.reform$score, bkgnd.phylo.reform$score)))
stats = append(stats, paste("\n\n"))
stats = append(stats, capture.output(wilcox.test(CL.phylo.reform$score, bkgnd.phylo.reform$score)))

writeLines(stats, "phyloP_CL.v.unbound_ks.test_wilcox.test_190408")




## get mean and median info for datasets
stats2= paste("Summary of CL nucleotide dataset")
stats2 = append(stats2, paste("\n"))
stats2 = append(stats2, capture.output(summary(CL.phylo.reform)))
stats2 = append(stats2, paste("\n\n\n\n\n\n"))
stats2= append(stats2, paste("Summary of noncrosslinked nucleotide dataset"))
stats2 = append(stats2, paste("\n"))
stats2 = append(stats2, capture.output(summary(bkgnd.phylo.reform)))

writeLines(stats2, "phyloP_CL.unbound_summary.stats_190408")


