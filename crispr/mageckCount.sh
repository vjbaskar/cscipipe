# Crispr analysis using mageck

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments
-f|--inputFile) [samples.mageck.txt]
-o|--organism) crispr organism specific library file eg. ~/GENOMES/CRISPR/*_mageck.csv
-0|--cmd0) trimLength [  Number to trim the 5 prime end. Or give it as AUTO ]
-1|--cmd1) "--test-run --unmapped-to-file unmapped.txt" [  Number to trim the 5 prime end. Or give it as AUTO ]

++ recommended: 
-F|--farm = flag
-m|--mem=2000   # becomes bigger with bam file size
-q|--queue=normal *If samples > 50 use basement*
-p|--procs=1


++ inputFile
<fastq.gz> <samplename> <condns>
fq1.gz fq1_sampleA ctrl
fq2.gz fq2_sampleB ctrl
fq3.gz fq3_sampleC test
fq4.gz fq4_sampleD test



EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile organism cmd0 mconf_mageck`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

trimLength=$cmd0
sgRNALib=`readlink -f $organism`

testMode=0
if [ ! -z ${cmd1} ]; then
	testMode=`echo $cmd1 | grep  "test" | wc -l`
fi


##### Get the file names
filenames=`awk ' { print $1 } ' $inputFile | xargs`
exptnames=`awk ' { print $2 } ' $inputFile | xargs | sed -e 's/ /,/g'`
conditions=`awk ' { print $3 } ' $inputFile | xargs | sed -e 's/ /,/g'`


## Count data
command="$mconf_mageck/mageck count -l $sgRNALib -n $title --sample-label \"$exptnames\" --fastq $filenames --trim-5 $trimLength --pdf-report $cmd1"


totalExpts=`cat $inputFile | wc -l`
if [ $testMode -eq 0 ]; then
	sms "Running in full mode"
	if [ $totalExpts -gt 10 ]; then
		warnsms "Looks like the total expts exceeds 10. Switching to long queue"
		export queue="long"
	fi
	if [ $totalExpts -gt 50 ]; then
		warnsms "Looks like the total expts exceeds 10. Switching to basement queue"
		export queue="basement"
	fi
	
fi
echo $command
# farm submission: commandfile var = commandFile

if [ ! -z "$farm" ]; then
	export command="$command"
    echo $command
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $command"
fi
