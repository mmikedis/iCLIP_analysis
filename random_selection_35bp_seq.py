## random_selection_35bp_seq.R
## 
## Maria Mikedis mikedis@wi.mit.edu	

## This script randomly selects 20 35bp sequences per fasta entry
## First input file is a fasta of 3' UTR exons (if a transcript's 3'UTR has multiple exons, the exons have already been merged)
## first line: >NM_001163314_utr3_0_0_chr1_54473000_r(-)
## second line: sequence, oriented in 5' to 3' direction for the associated strand




###run as:
	### python random_selection_21bp_seq.py input_file.fa output.fa

	

from __future__ import division
import sys
import os 
import re
import math
import random ### to use random number generate

###for division that needs to be done later



fasta = []


with open(sys.argv[1], 'r') as input:
#with open("mm10_GRCm38_MM257.tpm1_coding.utr3.exons_longest.utr3.per.gene_Dazl.bound.only_merged.fa", 'r') as input:	
	for line in input:
		mm = line.strip("\n")
		fasta.append(mm)

output = [] 


### randomly select 20 random 35bp sequences per fasta entry
i=0
while i<len(fasta):
	if fasta[i][0] != ">": ## if i divided by 2 does not produce a remainder of 0
		sequence = fasta[i]
		random_numbers = []
		for y in range(20): 	### randomly generate 20 numbers
			random_numbers.append(random.randint(0,len(sequence)-35))
		m=0
		random_numbers = list(set(random_numbers)) ### remove repeated integers from list
		while m<len(random_numbers):
			output.append((fasta[i-1])+"_random_seq"+str(m+1))	### output fasta header with random sequence numbered
			output.append(sequence[random_numbers[m]:(random_numbers[m]+35)])	### to get 35bp sequence starting at random number 
			m=m+1
	i=i+1

with open(sys.argv[2], 'w') as output_fasta:
	for line in output:
		new_line = "".join(str(v) for v in line)
		output_fasta.write(new_line+"\n")
