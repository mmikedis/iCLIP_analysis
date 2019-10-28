#### extract_3UTR_sequences.R
## 
## Maria Mikedis mikedis@wi.mit.edu	
## Run as: extract_3UTR_sequences.R ${output_filename}
## Extract the 3' UTR sequences from mm10 using UCSC table "refGene"
## Output 2 files
## File 1: ${output_filename}.fa of sequences
## File 2: ${output_filename}_length.txt containing lengths of 3'UTRs
	## format: Refseq_id 	length_3UTR
## Dependencies: GenomicFeatures, BSgenome.Mmusculus.UCSC.mm10


args <- commandArgs()
print(args)

filename <- args[6]  ### filename is passed as argument #6

library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)

## Create gene model of mm10 genome from UCSC table ncbiRefSeqCurated
txdb <- makeTxDbFromUCSC(genome="mm10", tablename="refGene")

## Load mm10 genome as DNAstring object
genome <- BSgenome.Mmusculus.UCSC.mm10

## Extract the transcript coordinates from this gene model
transcripts(txdb)

## Create a map between gene and transcript identifiers by outing both
txbygene <- transcriptsBy(txdb, "gene")

## Extract the 3' UTR coordinates; use.names option keeps the RefSeq ID with the coordinate
threeUTR <- threeUTRsByTranscript(txdb, use.names=TRUE)

## Extract the 3'UTRs from the BSgenome
threeUTR_seqs <- extractTranscriptSeqs(genome, threeUTR)

## Export sequences as fasta file
writeXStringSet(threeUTR_seqs, file=paste(filename, ".fa",sep = ""))

## Create dataframe containing transcript name and width of 3' UTR; write this dataframe to a file
df=data.frame(names(threeUTR_seqs), width(threeUTR_seqs)) 
colnames(df)=c("id", "3UTR_length")
write.table(df, file=paste(filename, "_length.txt",sep = ""), quote=F, sep="\t", row.names=F, col.names=T)

