# converting_mapped_reads_to_crosslinked_sites.py
# usage:
# python converting_mapped_reads_to_crosslinked_sites.py input.sam output.sam
# Maria Mikedis mikedis@wi.mit.edu

"""
# Take mapped CLIP data, input as a sam file
# Convert sam file information so that it only includes the crosslinked site, i.e., the nucleotide immediately 5' to the 5' start of the read
# Out sam file that contains crosslinked sites only
"""

from __future__ import division
import sys
import re
import os


with open(sys.argv[1], 'r') as input, open(sys.argv[2], 'w') as output:
	for line in input:
		field = line.split("\t")
		if field[1] == '0':  ### for reads mapping to forward strand, i.e., SAM bitwise flag (column 2) = 0
			field[3] = str(int(field[3]) - 1) ### the crosslinked position is the nucleotide immediately 5' to the 5' end of the read, i.e., position flag (column 4, representing leftmost mapping position) minus 1
			field[5] = "1M" ### changing read length to 1
			field[9] = field[9][0] ### modify query sequence so that it is only 1 nt (the nucleotide that comes after the crosslinked site)
			new_line = "\t".join(field[0:])
			new_line = str(new_line)
			output.write(new_line)
		if field[1] == '16': 	### for reads mapping to reverse strand, i.e., SAM bitwise flag (column 2) = 16
			read_length = len(field[9])  	### length of mapped read
			field[3] = str(int(field[3]) + read_length) #### the crosslinked position is the nucleotide immediately 5' to the 5' end of the read, i.e., position flag (column 4, representing leftmost mapping position) plus length of read [or position flag plus (length of read minus 1) plus +1]
			field[5] = "1M"
			field[9] = field[9][read_length-1] ### modify query sequence so that it is only 1 nt (the nucleotide that comes after the crosslinked site)
			new_line = "\t".join(field[0:])
			new_line = str(new_line)
			output.write(new_line)


