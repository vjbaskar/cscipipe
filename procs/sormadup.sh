Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Sorts, marks duplicates for sam, bam files

++ Mandatory args: 
-t|--title [project_Xmen]
-f|--inputFile [input file. see below]
-F|--farm [no input]
-m|--mem [10000]   # becomes bigger with bam file size
-p|--procs [ 1 ]
-P|--progArgs [optional for providing with other flags right out of the manual]
-o|--organism [ organism.fa ]
-0|--cmd0 [ input format: sam,bam ]
-1|--cmd1 [ output format: sam,bam,cram ]
-2|--cmd2 [ output sort order: coordinate, queryname ]

++ inputFile format 
a.bam
b.bam

++ Packages used
biobambam2/bamsormadup


EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title inputFile mem procs cmd0 cmd1 cmd2 organism`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

base_cmd="${mconf_biobambam2}/bamsormadup threads=${procs} inputformat=${cmd0} outputformat=${cmd1} SO=${cmd2} reference=${organism} "

# if [ ${dataType} == "nondirectional" ]; then
# 	infosms "Turning off directional alignment"
# 	base_cmd="${base_cmd} --non_directional "
# fi


export comments="sormadup"
commandFile="$title.${comments}.cmds"
rm -f $commandFile
echo -n "" > $commandFile
while read line
do
        l=($line)
        bam=${l[0]}
        b=`echo $bam | sed -e "s/${cmd0}$/sormadup.${cmd1}/g"`
        
        cmd="${base_cmd} M=${b}.metrics indexfilename=${b}.bai  < $bam > $b  "
        echo "${cmd}" >> $commandFile
        echo "${cmd}"
done < <(grep -v "^#" $inputFile)

# commands

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi
