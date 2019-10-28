#### obtain_3UTRseq_from_bed.R
## 
## Maria Mikedis mikedis@wi.mit.edu	


## run as: Rscript obtain_3UTRseq_from_bed.R ${input_bed} ${Refseq_ids} ${expand_regions} ${output_filename}
## ${input_bed} is a bed file containing regions within transcripts (i.e., crosslinked sites); chromosome field must be in UCSC style ("chr1")
	## format:
	## FDR must be in column #13
## ${Refseq_ids} is a list of transcripts (Refseq ids) from which to obtain sequence data for (filtered so that it contains the longest 3'UTR isoform per gene)
## ${expand_regions} is an integer by which the sequence originating from the bed file will be expanded right and left sides; expanded sequence will be the transcript sequence only (no genomic sequence)
## ${output_filename}
## uses mm10 GRCm38.90
## Dependencies: GenomicFeatures, BSgenome.Mmusculus.UCSC.mm10


##https://rdrr.io/bioc/GenomicFeatures/man/coordinate-mapping-methods.html


##############
### export relative position of crosslinked site within 3' UTR


args <- commandArgs()
print(args)

input <- args[6]  
input_ids <- args[7]
expand <- as.integer(args[8])
output <- args[9]  


library(BSgenome.Mmusculus.UCSC.mm10, lib.loc = "/home/mikedis/R/x86_64-pc-linux-gnu-library/3.5/")
library(GenomicRanges, lib.loc = "/home/mikedis/R/x86_64-pc-linux-gnu-library/3.5/")
library(ensembldb)
library(GenomicFeatures)
sessionInfo()
#lib.loc = "/home/mikedis/R/x86_64-pc-linux-gnu-library/3.5/"


## Load bedfile as genomic ranges (GRanges)
bedfile = read.table(input)
#bedfile = read.table("../background_transcriptome/mm10_GRCm38_MM257.tpm1_coding.utr3.exons_longest.utr3.per.gene.bed")
my.granges <- GRanges(seqnames = bedfile$V1, ranges = IRanges(start = (bedfile$V2), end = bedfile$V3, names = bedfile$V4), strand = bedfile$V6)

## Load list of transcript RefSeq Ids for which to get sequence information from
ids = read.table(input_ids)

## Create TxDb object from UCSC Genome Bioinformatics "RefGene" transcript table
supportedUCSCtables(genome="mm10")
txdb <- makeTxDbFromUCSC(genome="mm10", tablename="refGene")

## verify that style matches between TxDb and GRanges
seqlevelsStyle(txdb)
seqlevelsStyle(my.granges)

## Extract transcripts grouped by gene from TxDb object 
#transcriptsBytx = transcriptsBy(txdb, "gene")


## Extract the 3' UTR coordinates; use.names option keeps the RefSeq ID with the coordinate
threeUTR <- threeUTRsByTranscript(txdb, use.names=TRUE)

## Extract transcript coordinates for input GRanges
#txcoor=mapToTranscripts(my.granges, transcripts(txdb, use.names=TRUE), ignore.strand = FALSE)
txcoor=mapToTranscripts(my.granges, threeUTR, ignore.strand = FALSE)


## subset GRanges list so it contains my transcripts of interest
my.txcoor = txcoor[seqnames(txcoor) %in% as.vector(as.character(ids$V1))]



## Expand start and end of transcript coordinates
## trim function trims the coordinates that are out of range
txcoor_expanded = trim(my.txcoor + expand, use.names=TRUE)


## Convert transcript position info back to genomic position
gcoor = mapFromTranscripts(txcoor_expanded, threeUTR, ignore.strand = FALSE) 

names(gcoor)= paste(seqnames(gcoor), ":", start(ranges(gcoor)), "-", end(ranges(gcoor)), "(", strand(gcoor), ")", sep = "")


## edit Granges of gcoor so that UCSC format is maintained (txcoor_expanded doesn't follow UCSC format, so the formatting is disrupted when converting txcoor_expanded to gcoor)
## if sequence is on the positive strand, subtract 1 from end of range
## if sequence is on the negative strand, add 1 to the start of the range
i=1
len = length(strand(gcoor))
while (i<=len) 
{
	temp = (strand(gcoor)[i] == "+") 
	if (runValue(temp) == TRUE) 
	{
		end(ranges(gcoor))[i] = (end(ranges(gcoor))[i])-1
	}
	temp2 =  (strand(gcoor)[i] == "-") 
	if (runValue(temp2) == TRUE) 
	{
		start(ranges(gcoor))[i] = (start(ranges(gcoor))[i])+1
	}
	i=i+1
}



## Load mm10 genome as DNAstring object
genome <- BSgenome.Mmusculus.UCSC.mm10

## Extract transcript sequences
seqs = getSeq(Mmusculus, gcoor)

## Export sequences as fasta file
writeXStringSet(seqs, file=paste(output, ".fa" , sep = ""))
