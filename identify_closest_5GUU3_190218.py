#### identify_closest_5GUU3_190218.py
## 
## 2/18/2019
## Maria Mikedis mikedis@wi.mit.edu	
## run as: python identify_closest_5GUU3_190218.py ${input_fa} ${input_bed} ${output_filename}
## identifies the closest GUU to the crosslinked site for sequences plus and minus 10nt around crosslinked site

## ${input_fa} is a fasta 
	## first line: >chr7:116104561-116104584(+)
	## second line: GTGTTAAAATGTTTTGCAAAT
	## fasta file obtained from bed file with each line representing a crosslinked nucleotide
	## sequence is 21nt in length, with crosslinked site in position 11
	## sequences are in the correct 5' to 3' orientation based on strand

## ${input_bed} is a bed file containing crosslinked sites (1 nt per line)

## ${output_filename}
	## format: chr7 	start(CLnt-10) 	stop(CLnt-10) 	strand 	sequence 	start.of.GUU	stop.of.GUU 	orientation.of.GUU 	distance.GUU.from.CL.site 	-log10(pval)
	## for distance values, positive means GUU is 3' to the crosslinked site; negative means GUU is 5' to crosslinked site



from __future__ import division
import sys
import os 
import re
import math

###for division that needs to be done later


crosslinked = []

## load fasta file
with open(sys.argv[1], 'r') as crosslinked_sites:
	for line in crosslinked_sites:
		mm = line.strip("\n")
		split_line = mm.split("\t")
		crosslinked.append(split_line)


## load bed file
CLIP_scores = []
with open(sys.argv[2], 'r') as file:
	for line in file:
		mm = line.strip("\n")
		split_line = mm.split("\t")
		CLIP_scores.append(split_line)


## combine bed file info with fasta sequence info
input=[]
i=0
while i<len(crosslinked):
	if i%2==0:
		new_line = []
		read_line = str(crosslinked[i])
		chromosome = re.search(r'>(.*):', read_line)
		start = re.search(':(.*?)-', read_line)
		stop = re.search(r"-(.*?)\(", read_line)
		strand = re.search(r'\((.*?)\)', read_line)
		if chromosome:
			new_line.append(chromosome.group(1))
		if start:
			new_line.append(start.group(1))
		if stop:
			new_line.append(stop.group(1))
		if strand:
			new_line.append(strand.group(1))
		i=i+1
	else:
		new_line.append('\t'.join(crosslinked[i]))
		input.append(new_line)
		i=i+1


## identify closest GUU to each crosslinked site; if no GUU within plus or minus 10nt of the crosslinked site, return "n/a"
## check if crosslinked site (position 11) is a G; if not, move to right and left to find closest G
## all coordinates for GUU positioning are in BED format

output = []
i=0
while i<len(input):
	read_line = input[i]
	sequence=str(read_line[4])
	y=0	###variable for what position in sequence we are looking at for G
	if len(sequence)!= 21: 		### if the sequence is less than 21nt long, skip it and more to next sequence
		y=12
		i=i+1
	while y==0:
		if read_line[3]=="+":
			if sequence[(y+10):(y+13)] == str("GTT"):
				start = int(read_line[1])+int(10)
				end = int(read_line[1])+int(12+1)
				orientation = str("forward")
				distance=0
				read_line.append(start)
				read_line.append(end)
				read_line.append(orientation)
				read_line.append(distance)
				output.append(read_line)
				i=i+1 
				y=12
		if read_line[3]=="-":
			if sequence[(y+10):(y+13)] == str("GTT"):
				if read_line[3]=="+":		
					start = int(read_line[1])+int(8)
					end = int(read_line[1])+int(10+1)
				if read_line[3]=="-":		##in a sequence from negative strand, start value corresponds to last nt in sequence
					start = int(read_line[1])+int(10) ###position that matches 2nd T on negative strand
					end = int(read_line[1])+int(12+1) ###position that matches G on negative strand, plus one
				distance=0
				orientation = str("reverse")
				read_line.append(start)
				read_line.append(end)
				read_line.append(orientation)
				read_line.append(distance)
				output.append(read_line)
				i=i+1
				y=12
		y=y+1
	while y<11:
		if sequence[10-y] == str("G") and sequence[11-y]==str("T") and sequence[12-y]==str("T"):
			if read_line[3]=="+":		
				start = int(read_line[1])+int(10-y)
				end = int(read_line[1])+int(12-y+1)
			if read_line[3]=="-":		##in a sequence from negative strand, start value corresponds to last nt in sequence
				start = int(read_line[1])+int(8+y) ###position that matches 2nd T on negative strand
				end = int(read_line[1])+int(10+y+1) ###position that matches G on negative strand, plus one
			orientation = str("forward")
			if 11-y == 10 or 12-y == 10:
				distance = 0
			else: 
				distance = (2-y)
			read_line.append(start)
			read_line.append(end)
			read_line.append(orientation)
			read_line.append(distance)
			output.append(read_line)
			i=i+1 
			y=12
		elif (y+12)<21 and sequence[10+y] == str("G") and sequence[11+y]==str("T") and sequence[12+y]==str("T"):
			if read_line[3]=="+":		
				start = int(read_line[1])+int(10+y)
				end = int(read_line[1])+int(12+y+1)
			if read_line[3]=="-":		##in a sequence from negative strand, start value corresponds to last nt in sequence
				start = int(read_line[1])+int(8-y) ###position that matches 2nd T on negative strand
				end = int(read_line[1])+int(10-y+1) ###position that matches G on negative strand, plus one
			orientation = str("forward")
			distance=(y)
			read_line.append(start)
			read_line.append(end)
			read_line.append(orientation)
			read_line.append(distance)
			output.append(read_line)
			i=i+1
			y=12
		else:
			y=y+1
	if y==11:
		start = str("n/a")
		end = str("n/a")
		orientation = str("n/a")
		distance = str("n/a")
		read_line.append(start)
		read_line.append(end)
		read_line.append(orientation)
		read_line.append(distance)
		output.append(read_line)
		i=i+1


## for each GUU and its crosslinked site, add the -log base 10 of pvalue to end of line
i=0
while i<len(output):
	chromosome = output[i][0]
	crosslinked_start = int(output[i][1])+10 ### the output file contains the start site of the 21nt around crosslinked nucleotide; actual crosslinked nucleotide is start site + 11
	n=0
	for line in CLIP_scores:
		if line[0] == chromosome and line[1] == str(crosslinked_start):
			log = math.log(float(line[4]), 10)*-1
			output[i].append(log)
	i=i+1


## export output file
with open(sys.argv[3], 'w') as output_bed:
	for line in output:
		new_line = "\t".join(str(v) for v in line)
		output_bed.write(new_line+"\n")
