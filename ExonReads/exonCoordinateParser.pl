#!/usr/bin/perl

# Parses the coordinates of exons from the exonRef.txt file

open (FILE, 'exonRef_CCDS.txt');
LINE: while (<FILE>) {
chomp;
my @values;
($chromosome, $strand, $exonStart, $exonEnd) = split("\t");
next LINE if ($chromosome eq 'chrom');
# Form Arrays for exon start and end coordinates
@startCoords = split(',', $exonStart);
@endCoords = split(',', $exonEnd);

# Add all values in this line to recurring values
$siz = scalar(@startCoords);
# commenting out
#for ($z = 0; $z < $siz; $z++) {
# for ($y = $startCoords[$z]; $y <= $endCoords[$z]; $y++) {
# $temp = substr($chromosome,3);
#print "$chromosome $temp $strand $y \n";
#}
#}

# Instead of printing out individual coordinates, print out chr#:start:finish
for ($z = 0; $z < $siz; $z++) {
print "$chromosome:$startCoords[$z]:$endCoords[$z],$strand\n";
}
}
close (FILE);
exit;

