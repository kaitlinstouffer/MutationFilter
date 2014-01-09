#!/usr/bin/perl

# Parses the coordinates of exons from the exonRef.txt file

open (FILE, 'exonRefSeqFormatted_noDup.txt');
while (<FILE>) {
chomp;
($chromosome, $exonStart, $exonEndStrand) = split(":");
$chrNum = substr $chromosome, 3;
($exonEnd, $strand) = split(',', $exonEndStrand);

for ($z = int($exonStart); $z <= int($exonEnd); $z++) {
  print "$chrNum:$z\n";
}
}

close (FILE);
exit;
