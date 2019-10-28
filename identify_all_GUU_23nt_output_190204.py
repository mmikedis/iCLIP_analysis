#### identify_all_GUU_23nt_output_190204.py
### 2/4/19
### Maria Mikedis mikedis@wi.mit.edu  

### input is entire 3'UTR sequence in fasta format
### output is 5' GUU 3' sequences with a 23bp sequence total (3nt of GUU motif, plus and minus 10 nt)

###run as:
	### python identify_all_GUU_23nt_output_190204.py input_file.fa output.txt
	###output file format: fasta_header	sequence	start_of_GUU	stop_of_GUU	orientation	nt_from_CL_site 
		###for distance values, positive means GUU is 3' to the crosslinked site; negative means GUU is 5' to crosslinked site; for CLIP, expecting values to be negative, ie
	####positional output in bed file format

	
from __future__ import division
import sys
import os 
import re
import math
###for division that needs to be done later


fasta_data = []

with open(sys.argv[1], 'r') as fasta:
#with open("random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa", 'r') as fasta:
	for line in fasta:
		mm = line.strip("\n")
		split_line = mm.split("\t")
		fasta_data.append(split_line)



input=[]
i=0
while i<len(fasta_data):
	if i%2==0:
		new_line = []
		read_line = str(fasta_data[i])
		new_line.append(read_line)
		i=i+1
	else:
		new_line.append('\t'.join(fasta_data[i]))
		input.append(new_line)
		i=i+1


###find GUUs; output all GUUs, with 10nt on either side (G is nt 0, minus 10nt, plus 12nt)
###all coordinates for GUU positioning are in BED format
output = []	### for sequences with GUU
output2 = [] ### for sequences without any GUUs
i=0
while i<len(input):
	read_line = input[i]
	sequence=input[i][1]
	y=10	###variable for what position in sequence we are looking at for G; starting off 10 nt away from beginning of 3' UTR
	while y<len(sequence)-12:
		if sequence[y:(y+3)] == str("GTT"):	
			new_line = [] ### delete content of previous new_line list
			new_line.append(read_line[0]) ### add sequence header to new line
			new_line.append(sequence[(y-10):(y+13)]) ### add sequence to new line 
			#print(new_line)
			output.append(new_line)
		else: ## if no GUU
			output2.append(read_line) ### add line to output2; not writing this to file later; only using this for troubleshooting script
		y=y+1
	i=i+1

with open(sys.argv[2], 'w') as output_file:
	for line in output:
		new_line = "\t".join(str(v) for v in line)
		output_file.write(new_line+"\n")
