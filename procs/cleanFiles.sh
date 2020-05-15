# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Uses BWA aligner to align reads. 
Gives you sorted bam files.

++ Mandatory args: 
-t|--title [title]
-d|--dataType file extension [sam,cram,bam,fq.gz]
-P|--

++ inputFile format
2905STDY6551601 A1_DMSO1
2905STDY6551602 A2_DMSO2
2905STDY6551603 A3_DRUG1
2905STDY6551605 A4_DRUG2

++ Packages used
samtools, bwa

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



# farm submission: commandfile var = commandFile
# if [ ! -z "$farm" ]; then
# 	export fileOfCommands=$commandFile
# 	$mconf_installdir/bin/farmsub.sh
# else 
# 	warnsms "Not running as farm: Not recommended"
# 	warnsms "$mconf_bashexec $commandFile"
# fi
