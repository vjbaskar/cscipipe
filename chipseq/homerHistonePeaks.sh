# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments
-f|--inputFile

++ recommended: 
-F|--farm = flag
-m|--mem=2000   # becomes bigger with bam file size
-p|--procs=1

organism options: hs,mm

***** inputFile format
<tagDir test> <tagDir control>
homerTG_MOLM13_KAT7	homerTG_MOLM13_Input
homerTG_MV411_KAT7	homerTG_MV411_IgG
	
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
	command=""
    l=($line)
    treatment=${l[0]}
    control=${l[1]}
    t=`echo $treatment | sed -e 's/homerTG_//g'`
    c=`echo $control | sed -e 's/homerTG_//g'`
    $mconf_homer/findPeaks $treatment -o homerPeaks_${t}_.vs._${c} -style histone -i $control -fdr 0.05
done < $inputFile


# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi


