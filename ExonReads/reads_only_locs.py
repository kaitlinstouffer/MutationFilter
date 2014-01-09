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
# Assume line is in format of: chr#:location
# Don't need strand information; just tally up total reads
for line in exonCoords:
    wholeInput = line.split(':')
    position = int(wholeInput[1]) # assume coordinates are 1-based from vcf
    chromosome = wholeInput[0]# Get only part that has number if in form of chr#_contig

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
    
    region_string = chromosome + ":" + str(position) + ":" + str(position + 1) 
    iter = samfile.pileup(region=region_string,max_depth=20000);
    count = 0

    for x in iter:
     #   index = 0;
        # Ignore if outside exon bounds
	if (x.pos != (position-1)):
	    continue;
        count += 1

    # Check for whether some positions have no reads / coverage
    if (count == 0):
        print line

