#!/usr/bin/perl

# Parses the coordinates of exons from the input file and writes 
# the coordinates in the form of:
#	GENENAME
#	chr#	
#	start1:end1;start2:end2;start3:end3...
# Usage: perl ./exonCoordinateParser.pl fileNameIn
# Note: the same genes will be printed out with different versions of same
# chromosome

use strict;
use warnings;

open (FILE, $ARGV[0]);
LINE: while (<FILE>) {
	chomp;
	my @values;
	(my $name, my $chromosome, my $exonStart, my $exonEnd) = split("\t");
	# Skip Header Line
	next LINE if ($chromosome eq 'chrom');

	# Form Arrays for exon start and end coordinates
	my @startCoords = split(',', $exonStart);
	my @endCoords = split(',', $exonEnd);

	# Add all values in this line to recurring values
	my $siz = @startCoords;

	# Make Array of start:end coords
	my @coords;
	for (my $i = 0; $i < $siz; $i++) {
		push @coords, "$startCoords[$i]:$endCoords[$i]"
	}

	print "$name\n";
	print "$chromosome\n";
	print join(";", @coords) . "\n";

	# commenting out
	#for ($z = 0; $z < $siz; $z++) {
	# for ($y = $startCoords[$z]; $y <= $endCoords[$z]; $y++) {
	# $temp = substr($chromosome,3);
	#print "$chromosome $temp $strand $y \n";
	#}
	#}

	# Instead of printing out individual coordinates, print out chr#:start:finish
	#for ($z = 0; $z < $siz; $z++) {
	#	print "$chromosome:$startCoords[$z]:$endCoords[$z],$strand\n";
	#}
}
close (FILE);
exit;

