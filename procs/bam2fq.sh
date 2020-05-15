#!/bin/sh

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments
-f|--inputFile
-r|--paired = [0,1]

++ recommended: 
-F|--farm = flag
-m|--mem=20000   # becomes bigger with bam file size
-p|--procs=1

***** inputFile format
DSMX47.4g.sort.bam 
DSMX47.4h.sort.bam 
IHMX1.1c.sort.bam 
IHMX2.1b.sort.bam
EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile paired mconf_biobambam2`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi


bam2fastq() {
	INPUT=$1
	sampleID=`basename $INPUT .cram`
	sampleID=`basename $sampleID .bam`
	if [ $paired -eq 0 ]; then
	OUTPUT="$sampleID.fq"
	#tfile=`tempfile -d ./`
	command="$mconf_biobambam2/bamtofastq inputformat=bam S=$OUTPUT.gz filename=$INPUT gz=-1"
	echo "$command"
	fi

	if [ $paired -eq 1 ]; then
		OUTPUT0="${sampleID}.pair1.fq"
		OUTPUT1="${sampleID}.pair2.fq"
		command="$mconf_biobambam2/bamtofastq inputformat=bam F=$OUTPUT0.gz F2=$OUTPUT1.gz filename=$INPUT S=${sampleID}.single.fq.gz gz=-1"
		echo "$command"
	fi
}



commandFile="${title}.${comments}.cmds"
echo -n "" > $commandFile
while read line
do
	l=($line)
	id=${l[0]}
	name=`echo $id | sed -e 's/.cram//g' -e 's/.bam//g'`
	if [ "$paired" == 1 ]; then
		bam2fastq $id
	else 
		bam2fastq $id
	fi >>  $commandFile #$id.bam2fq.sh 
	#echo "sh $id.bam2fq.sh" >> $commandFile
done < $inputFile

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi

#farm..sub -q normal -n bam2fq -s $runName -m $MEM -p $THREADS -t 40 -a ${runName}.bam2fq.commands
