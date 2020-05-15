# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments
-f|--inputFile
-o|--organism 

++ recommended: 
-F|--farm = flag
-m|--mem=2000   # becomes bigger with bam file size
-p|--procs=1

organism options: hg38, hg19, mm10

***** inputFile format
sample1 sample1-1.bam,sample1-2.bam
sample2 sample2.bam
sample3 sample3-1.bam,sample3-2.bam
	
EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile organism`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi


commandFile="$title.$comments.cmds"
echo -n "" > $commandFile
while read line
do
	command=""
    l=($line)
    name=${l[0]}
    bams=`echo ${l[1]} | sed -e 's/,/ /g'`
    # convert bams to sams
    sams=""
    for b in $bams
    do
    	s=`echo $b | sed -e 's/.bam$/.sam/g'`
    	command="$command samtools view $b > $s ; "
    	sams="$sams $s"
    done
    
    command="$command $mconf_homer/makeTagDirectory homerTG_${name}/ $sams -format sam -genome ${organism} ; "
    command=" $command rm -f $sams " 
    echo $command >> $commandFile
    
done < $inputFile

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi
