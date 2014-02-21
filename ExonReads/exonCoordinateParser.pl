#!/usr/bin/perl

# Parses the coordinates of exons from the input file and writes 
# the coordinates in the form of:
#	GENENAME
#	chr#	
#	start1:end1;start2:end2;start3:end3...
# Usage: perl ./exonCoordinateParser.pl fileNameIn --coding 

# Note: --coding options specifies to consider only exon coordinates that are within
# coding regions of genes; genes with no coding regions will be eliminated
# Note: the same genes will be printed out with different versions of same
# chromosome

use strict;
use warnings;

my $CODING = "False";
if (@ARGV > 1) {
	if ($ARGV[1] eq "--coding") {
		$CODING = "True";
	}
}
		
my %genes;
open (FILE, $ARGV[0]);
LINE: while (<FILE>) {
	chomp;
	my @values;
	my ($name, $chromosome, $exonStart, $exonEnd, $cdsStart, $cdsEnd);
	if ($CODING eq "True") {
		($name, $chromosome, $exonStart, $exonEnd, $cdsStart, $cdsEnd) = split("\t");
		$cdsEnd =~ s/\r//g;
		next LINE if ($cdsStart eq $cdsEnd);
	}
	else {
		($name, $chromosome, $exonStart, $exonEnd) = split("\t");
	}

	# Skip Header Line
	next LINE if ($chromosome eq 'chrom');

	# Form Arrays for exon start and end coordinates
	my @startCoords = split(',', $exonStart);
	my @endCoords = split(',', $exonEnd);

	# Add all values in this line to recurring values
	my $siz = @startCoords;

	# Make Array of start:end coords
#	my @coords;
#	for (my $i = 0; $i < $siz; $i++) {
#		push @coords, "$startCoords[$i]:$endCoords[$i]"
#	}

	# Add Array of start:end coords to Gene if don't exist
	for (my $i = 0; $i < $siz; $i++) {
		if ($CODING eq "True") {
			if ($startCoords[$i] < $cdsStart) {
				$startCoords[$i] = $cdsStart;
				if ($endCoords[$i] < $cdsStart) {
					next;
				}
			}
			if ($endCoords[$i] > $cdsEnd) {
				$endCoords[$i] = $cdsEnd;
				if ($startCoords[$i] > $cdsEnd) {
					next;
				}
			}
		}
		my $entry = "$startCoords[$i]:$endCoords[$i]";
		if (not exists $genes{$name}{$chromosome}{$entry}) {
			$genes{$name}{$chromosome}{$entry} = 1;
		}
	}
}

foreach my $g (sort keys %genes) {
	foreach my $chr (sort keys %{ $genes{$g} }) {
		print "$g\n";
		print "$chr\n";
		print join(";", keys %{ $genes{$g}{$chr}}) . "\n";
	}
} 

#	print "$name\n";
#	print "$chromosome\n";
#	print join(";", @coords) . "\n";

#}
close (FILE);
exit;

