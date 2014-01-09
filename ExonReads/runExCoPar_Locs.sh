##Set PBS options
#PBS -N exonParserJob
#PBS -l walltime=8:00:00
cd $PBS_O_WORKDIR
python ./reads_only_locs.py 12.bam ./1.testLocs.txt >> ./12.testLocs.out.txt

