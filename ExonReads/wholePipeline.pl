#!/usr/bin/perl

while (<>) {
	($arg1, $arg2, $arg3, $arg4) = split(" ");
	my $ret = system("python ./bam_snp_stats.py A060034.bam $arg1 $arg2 $arg3 $arg4 >> ./exonLocsNotSequenced_new.txt");
}
	
