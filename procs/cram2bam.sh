# Converts cram file to bam file

Help() { 
echo -en "\033[31m"
cat << EOM
-------------------------------------
converts cram to bam (sorts by coord)
-------------------------------------

++ Mandatory args: 
-t|--title
-c|--comments
-f|--inputFile

++ recommended: 
-F|--farm = flag
-m|--mem=1000   # becomes bigger with bam file size
-p|--procs=1

***** inputFile format
input1.cram
input2.cram

**** Packages used
samtools

EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

commandFile="$title.$comments.cmds"
echo -n "" > $commandFile
while read line
do
	a=($line)
	cramfile=${a[0]}
	f=`basename $cramfile .cram`
	#echo "samtools view -b $cramfile -o - | samtools sort - -o $f.bam -O BAM" >> $commandFile
	echo "samtools view -b $cramfile -o - | bamsort inputformat=bam outputformat=bam markduplicates=1 index=1 indexfilename=$f.bam.bai > $f.bam" >> $commandFile
done < $inputFile

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
 	export fileOfCommands=$commandFile
 	$mconf_installdir/bin/farmsub.sh
 else 
 	warnsms "Not running as farm: Not recommended"
 	warnsms "$mconf_bashexec $commandFile"
 fi
