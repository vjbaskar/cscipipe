Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Generates input files for methylKit from Bismark bam files

++ Mandatory args: 
-t|--title [project_Xmen]
-f|--inputFile [input file. see below]
-F|--farm [no input]
-m|--mem [10000]   # becomes bigger with bam file size
-o|--organism [ mm10,hg38 Just a name to identify. Not used ]

++ inputFile format (bam files should be sorted and indexed)
a.sormadup.bam
b.sormadup.bam

++ Packages used
R methylKit


EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title inputFile farm mem organism`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi



base_cmd="${mconf_R} --vanilla "

# if [ ${dataType} == "nondirectional" ]; then
# 	infosms "Turning off directional alignment"
# 	base_cmd="${base_cmd} --non_directional "
# fi


export comments="methylproc"
commandFile="$title.${comments}.cmds"
rm -f $commandFile
echo -n "" > $commandFile
while read line
do
        l=($line)
        bam=${l[0]}
        #b=`echo $bam | sed -e "s/${cmd0}$/sormadup.${cmd1}/g"`
        
        cmd="${base_cmd} --args $organism $bam < ${mconf_installdir}/bin/methylseq/methylKitPrelim.R"
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
