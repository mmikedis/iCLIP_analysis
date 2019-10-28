#### identify_all_tri.nt_motifs_fa_input_only_180409.py
## 
## 2/18/2019
## Maria Mikedis mikedis@wi.mit.edu	
## run as: python identify_all_tri.nt_motifs_fa_input_only_180409.py ${tri-nucleotide motif} ${input_fa} ${output_filename}
## outputs all 5' NNN 3' sequences relative to the center nucleotide in a 21nt sequence
## (if the 21bp sequences are from CLIP data, the center nucleotide is the crosslinked site)


## ${input_fa} is a fasta from UCSC genome browser; 
## first line: >header_any_format
## second line: GTGTTAAAATGTTTTGCAAAT
## sequence is 21nt in length, with crosslinked site in position 11

## ${output_filename}
	## format: fasta.header	sequence	start.of.NNN	stop.of.NNN	orientation	distance.from.CL.site
	## starts and stops are within each individual 21nt sequence; range is fully closed
	## positive distance means NNN is 3' to the crosslinked site 
	## negative distance means NNN is 5' to crosslinked site


from __future__ import division
import sys
import os 
import re
import math

## load tri-nucleotide motif
motif=sys.argv[1]

#check that length of motif is 3 nucleotides; otherwise return error
if len(motif)!=3:
	sys.exit("Error! Motif must be 3 nucleotides long")

#convert Us to Ts because tfasta file contains Ts
if motif[0]=="U":
	motif = 'T' + motif[1:]

if motif[1]=="U":
	motif = motif[0:1] +'T' + motif[2:]

if motif[2]=="U":
	motif = motif[0:2] + 'T'

crosslinked = []

## load fasta file
with open(sys.argv[2], 'r') as crosslinked_sites:
#with open("utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line.fa", "r") as crosslinked_sites:
	for line in crosslinked_sites:
		mm = line.strip("\n")
		split_line = mm.split("\t")
		crosslinked.append(split_line)


## reorganize fasta file so that header and sequence are a single line
input=[]
i=0
while i<len(crosslinked):
	if i%2==0:
		new_line = []
		read_line = str(crosslinked[i])
		new_line.append(read_line)
		i=i+1
	else:
		new_line.append('\t'.join(crosslinked[i]))
		input.append(new_line)
		i=i+1


## find all NNNs in each sequence and output position relative to the center nucleotide of each sequence
## check the center of the sequence first and then move to the right and left of center
output = []
i=0
while i<len(input):
	no_motif_found = True
	read_line = input[i]
	sequence=str(read_line[1])
	if len(sequence)==35: ##only proceed if sequence length is 35nt
		y=0	## variable for what position in sequence we are looking at for the first nucleotide in tri-nucleotide sequence NNN
		while y==0:
			if sequence[y+17] == str(motif[0]) and sequence[y+18]==str(motif[1]) and sequence[y+19]==str(motif[2]):	
				new_line = [] ## delete content of previous new_line list
				new_line.append(read_line[0]) ### add sequence header to new line
				new_line.append(read_line[1]) ### add sequence to new line 
				start = int(17)
				end = int(19)
				distance=0
				orientation = str("forward")
				new_line.append(start)
				new_line.append(end)
				new_line.append(orientation)
				new_line.append(distance)
				output.append(new_line)
				no_motif_found = False
			y=y+1
		while y<18:
			if sequence[17-y] == str(motif[0]) and sequence[18-y]==str(motif[1]) and sequence[19-y]==str(motif[2]):	
				new_line = [] ## delete content of previous new_line list
				new_line.append(read_line[0]) ## add sequence header to new line
				new_line.append(read_line[1]) ## add sequence to new line 
				start = int(17-y)
				end = int(19-y)
				orientation = str("forward")
				distance = (y*-1)
				new_line.append(start)
				new_line.append(end)
				new_line.append(orientation)
				new_line.append(distance)
				output.append(new_line)
				no_motif_found = False
			if (y+19)<35 and sequence[17+y] == str(motif[0]) and sequence[18+y]==str(motif[1]) and sequence[19+y]==str(motif[2]):		
				new_line = [] ## delete content of previous new_line list
				new_line.append(read_line[0]) ## add sequence header to new line
				new_line.append(read_line[1]) ## add sequence to new line 
				start = int(17+y)
				end = int(19+y)
				orientation = str("forward")
				distance=(y)
				new_line.append(start)
				new_line.append(end)
				new_line.append(orientation)
				new_line.append(distance)
				output.append(new_line)
				no_motif_found = False
			y=y+1
		if y==18 and no_motif_found == True:
			new_line = [] ## delete content of previous new_line list
			new_line.append(read_line[0]) ## add sequence header to new line
			new_line.append(read_line[1]) ## add sequence to new line 
			start = str("n/a")
			end = str("n/a")
			orientation = str("n/a")
			distance = str("n/a")
			new_line.append(start)
			new_line.append(end)
			new_line.append(orientation)
			new_line.append(distance)
			output.append(new_line)
	i=i+1



with open(sys.argv[3], 'w') as output_bed:
	for line in output:
		new_line = "\t".join(str(v) for v in line)
		output_bed.write(new_line+"\n")
