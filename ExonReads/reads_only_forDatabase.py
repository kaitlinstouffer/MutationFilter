#!/usr/bin/python

# Usage: python ./reads_only_forDatabase.py bamFileName geneFileName > outputFile.txt
# Output File in form of:
#	GENE	AVG_COVERAGE	%_BASES_ZERO	%_BASES_LESS_20	NUM_EXONS	NUM_BASES	Exon1	Exon2	Exon3...
# GENE = genename
# AVG_COVERAGE = sum(avg_coverage per exon) / # exons
# %_BASES_ZERO = sum(num_bases < 0 per exon) / num_bases total
# %_BASES_LESS_20 = sum(num_bases < 20 per exon)
# NUM_EXONS = # exons per gene
# NUM_BASES = sum(num bases per exon)
# ExonX = avg_coverage for exon; num bases = 0; num bases < 20; num bases
# 	  avg_coverage for exon = sum (coverage per base) / num bases
#
# Note: exons may overlap so bases may be counted more than once in num_bases or in percentages

import sys
import getopt
import pysam

bamfile = sys.argv[1]
samfile = pysam.Samfile( bamfile, 'rb' )

references = samfile.references
chromosomeRef = True

if (not references[0].startswith('c')):
	chromosomeRef = False

# print Header
print "GENE\tAVG_COVERAGE\t%_BASES_ZERO\t%_BASES_LESS_20\tNUM_EXONS\tNUM_BASES\tEXONS\n"

exonCoords = open(sys.argv[2], 'r')
line = exonCoords.readline()
while (line):
	# Find Gene (Assume file starts with gene name)
	geneDict = dict(name=line.strip(),avg_coverage=0.0, perBaseZero=0.0, perBaseLess20=0.0,num_exons=0, num_bases=0)
	# Find chromosome
	chromosome = ""
	line = exonCoords.readline()
	chromosome = line.strip()
	chromosome_chr = chromosome
	
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

	# Get rid of non standard chromosomes
	if (len(chromosome.split('_')) >= 3):
		line = exonCoords.readline()
		line = exonCoords.readline()
		continue;

	if (not chromosomeRef):
		chromosome = chromosome[3:]
	
	line = exonCoords.readline()
	
	# Parse ExonCoordinates
	line = line.strip()
	exonDict = dict() #keep list of coordinates per exon
	exons = line.split(";")
	exonCount = 0
	for x in exons:
		startEnd = x.split(":")
		temp = []
		for z in range(int(startEnd[0]), (int(startEnd[1]) + 1)):
			temp.append(z)
		# Print Exon coordinates as given in reference so can compare between files
		exonName = "Exon" + str(exonCount) + ": " + chromosome_chr + ":" + startEnd[0] + "-" + startEnd[1]
		exonDict[exonName] = temp
		exonCount += 1
	geneDict['num_exons'] = exonCount

	for key in exonDict:
		num_bases = len(exonDict[key])
		avg_cov_base = 0.0
		num_bases_zero = 0
		num_bases_l20 = 0
		for position in exonDict[key]:
			region_string = chromosome + ":" + str(position) + ":" + str(position + 1) 
    			iter = samfile.pileup(region=region_string,max_depth=20000);
    			count = 0

    			for x in iter:
             			# Ignore if outside exon bounds
				if (x.pos == (position-1)):
	    				count = x.n
					break
			if (count == 0):
				num_bases_zero += 1
			if (count < 20):
				num_bases_l20 += 1
			avg_cov_base += count
		avg_cov_base = (avg_cov_base / num_bases)
		# put information in geneDict
		geneDict[key] = str(avg_cov_base) + ";" + str(num_bases_zero) + ";" + str(num_bases_l20) + ";" + str(num_bases)
		geneDict['avg_coverage'] += avg_cov_base
		geneDict['perBaseZero'] += num_bases_zero
		geneDict['perBaseLess20'] += num_bases_l20
		geneDict['num_bases'] += num_bases

	geneDict['avg_coverage'] = geneDict['avg_coverage'] / geneDict['num_exons']
	geneDict['perBaseZero'] = 100 * (geneDict['perBaseZero'] / geneDict['num_bases'])
	geneDict['perBaseLess20'] = 100 * (geneDict['perBaseLess20'] / geneDict['num_bases'])

	# print gene information
	print geneDict['name'] + "\t" + str(geneDict['avg_coverage']) + "\t" + str(geneDict['perBaseZero']) + "\t" + str(geneDict['perBaseLess20']) + "\t" + str(geneDict['num_exons']) + "\t" + str(geneDict['num_bases'])
	for key in geneDict:
		if (key.startswith("Exon")):
			print "\t" + key + "\t" + geneDict[key]
	print "\n"
	line = exonCoords.readline()		
