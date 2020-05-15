# Script: MACS2 data

# Main config file
#source ~/.$USER/cscipipe.conf
#! $mconf_bashexec

# get the exec dir
# basedir=`readlink -f $0 | xargs dirname`
#source $mconf_installdir/src/general_args.sh
#source $mconf_installdir/src/basic_functions.sh

#if [ ! -z "$conf" ]; then
#	source $conf # config file in the local directory can override the default ones
#fi


if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: 
	title
	comment
	mem
	proc
	inputFile
	
	++ Recommended to run in farm:
	-F|farm
	-m|memory 500
	-p|procs 1
	
	submits as array since $inputFile is a text file
	
	***** inputFile format
	bamfileA bamfile1.bam,bamfile2.bam,bamfile3.bam
	bamfileB bamfile4.bam,bamfile2.bam

	# Output is SO= coordinate and indexed
	"
	exit 2
fi

commandFile="$title.bammerge.cmds"
echo -n "" > $commandFile
while read line
do
    l=($line)
    name=${l[0]}
    files=${l[1]}
    fileList=`echo $files | sed -e 's/,/ I=/g'`
    fileList="I=$fileList"
	cmd="$mconf_biobambam2/bammerge $fileList SO=coordinate index=1 indexfilename=${name}.bam.bai > ${name}.bam"
	echo $cmd >> $commandFile
    #command="samtools merge - $fileList | samtools sort - -o $name.sort.bam ; samtools index $name.sort.bam"
    #farm..sub -n bammerge -s $name -p $threads -q long -m $mem -c "$command"
done < $inputFile

if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	export concurrentJobs=20
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	$mconf_bashexec $commandFile
fi


#bammerge I=MOLM13_KAT7_R1.bam I=MOLM13_KAT7_R2.bam SO=coordinate index=1 indexfilename=MOLM13_KAT7.bam.bai > MOLM13_KAT7.bam
