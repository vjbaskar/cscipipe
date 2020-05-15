# Converts cram file to bam file

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

organism options: hs,mm

***** inputFile format
input1.cram
input2.cram

**** Packages used
samtools
biobambam2

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
	echo "samtools view $cramfile -o - -b | samtools sort - -o - -O BAM | bamtofastq F=${f}.pair1.fq.gz F2=${f}.pair2.fq.gz S=${f}.fq.gz gz=1" >> $commandFile
done < $inputFile

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
 	export fileOfCommands=$commandFile
 	$mconf_installdir/bin/farmsub.sh
 else 
 	warnsms "Not running as farm: Not recommended"
 	warnsms "$mconf_bashexec $commandFile"
 fi