##Set PBS options
#PBS -N compareReads
#PBS -l walltime=8:00:00
cd $PBS_O_WORKDIR
if [ "x" == "x$FILE_IN" ] ; then
	echo "Variable is not set"
else
	sort ./$FILE_IN.txt -o ./$FILE_IN.sort.txt
	comm -1 -2 ./$FILE_IN.sort.txt ./exonRefSeq_IndividualLocs_sortUniq.txt > ./$FILE_IN.sort.intersect.txt
fi

