##Set PBS options
#PBS -N exonParserJob
#PBS -l walltime=8:00:00
cd $PBS_O_WORKDIR
##cat ./exonRefFormatted.txt | perl ./wholePipeline.pl
##perl ./exonCoordinateParser.pl | awk '!x[$0]++' > ./exonRefFormatted_ccds.txt
python ./reads_only.py 1.bam ./exonRefFormatted_ccds.txt >> ./exonLocsNotSequenced_ccds_1.txt

