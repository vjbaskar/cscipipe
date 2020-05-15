# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments
-f|--inputFile = samples file
-0|--cmd0 = comparisons file
-1|--cmd1  = input count file from mageckCount command (normalised - count_normalized.txt)

++ recommended: 
-F|--farm = flag
-m|--mem=1000   # becomes bigger with bam file size
-p|--procs=1
-q|--queue=small

***** inputFile = samples file format
<fastq.gz> <samplename> <condns>
fq1.gz fq1_sampleA ctrl
fq2.gz fq2_sampleB ctrl
fq3.gz fq3_sampleC test
fq4.gz fq4_sampleD test

***** cmd0 = comparisons file format
control1 cond2
control1 cond3
control2 cond2

	
EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile cmd0`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

countFile=$cmd1
comparisonsFile=$cmd0

getSamples() {
	condn=$1
	file=$2

	samples=`awk -v c=$condn ' { if($3==c) print $2 } ' $file | xargs | tr " " ","`
	echo $samples
}

commandFile="$title.$comments.cmds"
echo -n "" > $commandFile
while read line
do
	l=($line)
	expt1=${l[0]}
	expt2=${l[1]}
	samples1=`getSamples $expt1 $inputFile`
	samples2=`getSamples $expt2 $inputFile`
	
	n="mageck_${expt1}..vs..${expt2}"
	echo "mageck test -k $countFile -t $samples2 -c $samples1 -n $n --norm-method \"none\" --gene-test-fdr-threshold 0.20 --keep-tmp ; mkdir ${n} ; mv $n.* ${n} " >> $commandFile
done < $comparisonsFile

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi
