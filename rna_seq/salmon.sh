# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments
-f|--inputFile
-o|--organism : Look at $mconf_salmonindex
-0|--cmd0 : Library type [Use ISR for Ill-PE] - See below
-r|--paired: 0,1

++ recommended: 
-F|--farm = flag
-m|--mem=5000   # becomes bigger with bam file size
-p|--procs=1
-q|--queue=small


***** inputFile format
file1 file1.pair1.fq.gz file1.pair2.fq.gz
file2 file2.pair1.fq.gz file2.pair2.fq.gz


***** Library types
http://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype

--------------------------------------
TopHat	        Salmon (and Sailfish)
--------------------------------------
 	           	Paired-end	Single-end
-fr-unstranded		IU		U
-fr-firststrand		ISR		SR
-fr-secondstrand	ISF		SF


	
EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile organism paired cmd0`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

libType=$cmd0

#command="salmon quant -l $libType -i $organism" 

commandFile="$title.$comments.cmds"
echo -n "" > $commandFile
while read line
do
	command="$mconf_salmon/salmon quant --no-version-check -l $libType -i $organism" 
	l=($line)
	n=${l[0]}
	if [ "$paired" -eq 1 ]; then
		pair1=`echo ${l[1]} | sed -e 's/,/ /g'`
		pair2=`echo ${l[2]} | sed -e 's/,/ /g'`
		command="$command -1 $pair1 -2 $pair2"
	fi
	if [ "$paired" -eq 0 ]; then
		fq=`echo ${l[1]} | sed -e 's/,/ /g'`
		command="$command -r $fq"
	fi
	command="$command -o $n"
	echo $command >> $commandFile
done < $inputFile


	
 #-1 TCGA-AB-2806.pair1.fq.gz -2 TCGA-AB-2806.pair2.fq.gz -o TCGA-AB-2806_salmon




# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi
