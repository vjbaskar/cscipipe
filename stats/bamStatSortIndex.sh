Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title [title of run]
-f|--inputFile [inputfile]
-F|--farm [no args]
-m|--mem [20000]   # becomes bigger with bam file size
-p|--procs [1]

***** inputFile format
file1.bam
file2.bam

**** Packages used
biobambam2/bamsormadup

EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title inputFile farm mem procs`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi


cmd() {
inputBam=$1
OUTPUT_PREFIX=`basename $inputBam .bam`
echo "$mconf_biobambam2/bamsormadup SO=coordinate M=temp.metrics indexfilename=${OUTPUT_PREFIX}.sormadup.bam.bai threads=${procs} < ${inputBam} > ${OUTPUT_PREFIX}.sormadup.bam ; samtools flagstat ${OUTPUT_PREFIX}.sormadup.bam > ${OUTPUT_PREFIX}.sormadup.bamstats"
}

nameFile="$title.bamsortindex.cmds"
echo -n "" > $nameFile
while read line
do
		l=($line)
		id=${l[0]}
		cmd $id
done < $inputFile >>  $nameFile

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export fileOfCommands=$nameFile
	export comments="bamsortindex"
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi

