#### identify_all_GUU.UorA_fa_input_only_180409.py
### 4/9/18
### Maria Mikedis mikedis@wi.mit.edu  

### this script outputs all 5' GUU(U/A) 3' sequences with a 35bp sequence
### the output is the position of the G of a 5' GUU(U/A) 3' sequence relative to the middle nucleotide in a 35 bp sequence (if the 35bp sequences are from CLIP data, the middle nucleotide is the crosslinked site)

###input file is a fasta from UCSC genome browser; 
#####first line: >mm10_ct_iCLIPundiffgonia_3643_Lypla1 range=chr1:4845113-4845133 5'pad=10 3'pad=10 strand=+ repeatMasking=none
#####second line: GTGTTAAAATGTTTTGCAAAT
#####fasta file obtained from bed file with each line representing a crosslinked nucleotide
###sequence is 21nt in length, with crosslinked site in position 11


###run as:
	### python identify_closest_GUU.py input_file.fa output_with_GUU_position.bed
	###output file format: fasta_header	sequence	start_of_GUU	stop_of_GUU	orientation	nt_from_CL_site 
		###for distance values, positive means GUUG is 3' to the crosslinked site; negative means GUUG is 5' to crosslinked site; for CLIP, expecting values to be negative, ie
	####positional output in bed file format

	
from __future__ import division
import sys
import os 
import re
import math

###for division that needs to be done later


crosslinked = []

with open(sys.argv[1], 'r') as crosslinked_sites:
#with open("random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa", 'r') as crosslinked_sites:
	for line in crosslinked_sites:
		mm = line.strip("\n")
		split_line = mm.split("\t")
		crosslinked.append(split_line)



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


###check if crosslinked site (position 17) is a G; if not, move to right and left to find closest G
###all coordinates for GUU positioning are in BED format
output = []
i=0
while i<len(input):
	no_motif_found = True
	read_line = input[i]
	sequence=str(read_line[1])
	y=0	###variable for what position in sequence we are looking at for G
	while y==0:
		if len(sequence)>20 and sequence[y+17] == str("G") and sequence[y+18]==str("T") and sequence[y+19]==str("T") and (sequence[y+20]==str("T") or sequence[y+20]==str("A")):	
			new_line = [] ### delete content of previous new_line list
			new_line.append(read_line[0]) ### add sequence header to new line
			new_line.append(read_line[1]) ### add sequence to new line 
			start = int(17)
			end = int(20)
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
		if len(sequence)>(20-y) and sequence[17-y] == str("G") and sequence[18-y]==str("T") and sequence[19-y]==str("T") and (sequence[20-y]==str("T") or sequence[20-y]==str("A")):	
			new_line = [] ### delete content of previous new_line list
			new_line.append(read_line[0]) ### add sequence header to new line
			new_line.append(read_line[1]) ### add sequence to new line 
			start = int(17-y)
			end = int(20-y)
			orientation = str("forward")
			distance = (y*-1)
			new_line.append(start)
			new_line.append(end)
			new_line.append(orientation)
			new_line.append(distance)
			output.append(new_line)
			no_motif_found = False		
		if (y+20)<35 and len(sequence)>(20+y) and sequence[17+y] == str("G") and sequence[18+y]==str("T") and sequence[19+y]==str("T") and (sequence[y+20]==str("T") or sequence[y+20]==str("A")):		
			new_line = [] ### delete content of previous new_line list
			new_line.append(read_line[0]) ### add sequence header to new line
			new_line.append(read_line[1]) ### add sequence to new line 
			start = int(17+y)
			end = int(20+y)
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
		new_line = [] ### delete content of previous new_line list
		new_line.append(read_line[0]) ### add sequence header to new line
		new_line.append(read_line[1]) ### add sequence to new line 
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



with open(sys.argv[2], 'w') as output_bed:
	for line in output:
		new_line = "\t".join(str(v) for v in line)
		output_bed.write(new_line+"\n")
