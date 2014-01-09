#!/usr/bin/python
import sys
import getopt
import pprint
pp = pprint.PrettyPrinter(indent=4)

import pysam

bamfile = sys.argv[1]
samfile = pysam.Samfile( bamfile, "rb" )

# Check for format of reference (chr# or #); default is chr#
references = samfile.references
chromosomeRef = True

if (not references[0].startswith('c')):
	chromosomeRef = False
	

MIN_MAPQ         =  25
MIN_BASEQ        =  25


# Set appropriate variables based on assumed format to script of:
# chr# # strand position
# chromosome = '2';
exonCoords = open(sys.argv[2], 'r')
# Assume line is in format of: chr#:start:end,strand
# Don't need strand information; just tally up total reads
for line in exonCoords:
    wholeInput = line.split(',')
	    # strand_in = wholeInput[1]
    positionInfo = wholeInput[0].split(':')
    positionStart = str(int(positionInfo[1]) + 1)
    chromosome = positionInfo[0]# Get only part that has number if in form of chr#_contig

    # Deal with mitochondrial dna (put into format of MT)
    # Note that files that have chr# references have chrM, not chrMT (as of 29/11/2013)
    if (chromosome == 'chrM' and not chromosomeRef):
	chromosome = 'chrMT'

    indOfGl = chromosome.find('gl')
    if (indOfGl >= 0 and not chromosomeRef):
        indOfGl = indOfGl + 2
	indOfRand = chromosome.find('random')
	if (indOfRand >= 0):
	    indOfRand = indOfRand - 1
	    chromosome = 'chrGL' + chromosome[indOfGl:indOfRand] + '.1'
	else:
            chromosome = 'chrGL' + chromosome[indOfGl:] + '.1'

#    indOfCtg = chromosome.find('ctg')
    if (len(chromosome.split('_')) >= 3):
        continue;   
    
    if (not chromosomeRef):
	    chromosome = chromosome[3:]
    
    region_string = chromosome + ":" + positionStart + ":" + positionInfo[2]
    iter = samfile.pileup(region=region_string,max_depth=20000);
    counts = dict()
    for i in range(int(positionStart),int(positionInfo[2])+1):
        counts[i] = 0

    for x in iter:
     #   index = 0;
        # Ignore if outside exon bounds
	if (x.pos < int(positionStart) or x.pos > int(positionInfo[2])):
	    continue;
        counts[x.pos] += x.n

    # Check for whether some positions have no reads / coverage
    for x in counts.keys():
        if (counts[x] == 0):
	    print region_string + "\n"
	    break;


	




            

