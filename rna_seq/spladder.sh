Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Uses BWA aligner to align reads. 
Gives you sorted bam files.

++ Mandatory args: 
-t|--title [project_Xmen]
-c|--comments [spladder]
-f|--inputFile [input file. see below]
-g|--gtf [ gtf file : /lustre/scratch119/casm/team163gv/SHARED/reference/spladder/ ]
-F|--farm [no input]
-m|--mem [1000]   # becomes bigger with bam file size
-p|--procs [1]



++ inputFile format
f1.bam
f2.bam


++ Packages used
spladder.py


EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile farm mem procs gtf`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi



commandFile="$title.${comments}.cmds"
cmd="${mconf_spladder}/spladder.py -a ${gtf}
rm -f $commandFile
echo -n "" > $commandFile
while read line
do
        l=($line)
        id=${l[0]}
        name=${l[1]}

        if [ $paired -eq 1 ]; then
            echo "Running in paired end mode : $id"
            peAlign $id $name > $id.$name.bwa.sh
        else
            echo "Running in single end mode"
            seAlign $id $name > $id.$name.bwa.sh
        fi
        #sub_ncores normal bwa 20000 ${procs} $name "sh $id.$name.bwa.sh"
        echo "${mconf_bashexec} $id.$name.bwa.sh" >> $commandFile

done < <(grep -v "^#" $inputFile)

# commands

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	#$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi
