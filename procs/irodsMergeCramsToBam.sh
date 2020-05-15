# irods sample folders contain cram files. This program merges them into one.


Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments
-f|--inputFile

++ recommended: 
-F|--farm = flag
-m|--mem=1000   # becomes bigger with bam file size
-p|--procs=1


***** inputFile format
4010STDY7036819	C_3143NTA
4010STDY7036820	C_3143NTB
4010STDY7036821	C_3143NTC
4010STDY7036822	C_3143CISWTA
EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	Help
	warnsms "mandatory check failed"
	errorsms "Diagnostics:
	flags/parameters not defined: $mandatory_fails"
	
fi

#$mconf_Rscript $mconf_installdir/bin/createProject.R -s $namefile -c SAMPLEID,SAMPLENAME | sed -e '1d' > samples.mergecram.txt
#namefile="samples.mergecram.txt"


commandFile="$title.mergecram.cmds"
echo -n "" > $commandFile
while read line
	do
	l=($line)
	sampleID=${l[0]}
	sampleName=${l[1]}
	echo "Merging CRAM files from $sampleID to $sampleName"


	ls $sampleID/finalfiles/*.cram > ${sampleID}.tomerge.txt
	c=`cat ${sampleID}.tomerge.txt | wc -l`
	if [ "$c" -lt 1 ]; then
		echo "Error: No files found to merge"
		#exit 0
	fi

	if [ "$c" -eq 1 ]; then
		echo "Found only one file. No merging required"
		file=`cat ${sampleID}.tomerge.txt`
		#cp $file $sampleName.cram
		command="$mconf_samtools view -@ ${procs} -b -o $sampleName.bam $file"
		echo $command >> $title.mergecram.cmds
		#exit 0
	fi

	if [ "$c" -gt 1 ]; then
		command="$mconf_samtools merge -f -c -p -@ ${procs} -b ${sampleID}.tomerge.txt -O BAM  $sampleName.bam"
		echo $command >> $title.mergecram.cmds
	fi

done < $inputFile


#farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	
	export fileOfCommands=$commandFile
	#export concurrentJobs=20
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi

# cjobs=`cat $namefile.commands | wc -l`
# farm..sub -q normal -n mergecram -m $mem -p $threads -t 20  -s ${namefile} -a $namefile.commands
