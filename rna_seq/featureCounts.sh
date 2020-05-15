# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
FeatureCounts for generating count table.

++ Mandatory args: 
-t|--title [project_Xmen]
-c|--comments [fc]
-f|--inputFile [input file. see below]
-F|--farm [no input]
-m|--mem [20000]   # becomes bigger with bam file size
-p|--procs [1]
-r|--paired [0/1]
-o|--organism [org.fa file]
-g|--gtf [gene definition file]
-0|--cmd0 [strand: 0 = unstranded, 1 = illumina,fr-firststrand, 2 = solid,reverse]
-P|--progArgs [Additional args: " -Q 30 --ignoreDup -J -B "] See featureCounts for details.
-O|--outpt [counts.txt]

++ inputFile format
file1.bam
file2.bam


++ Packages used
subreads, featureCount

EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile paired organism gtf`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

if [ "$paired" -eq 1 ]; then
	pair_flag=" -p"
else
	pair_flag=""
fi

bams=`cat $inputFile | xargs`
total_bams=`cat $inputFile | xargs -n 1 | wc -l`

echo " total = $total_bams"
if [ ! -z "$outpt" ]; then
	warnsms "** Output file exists"
	outpt_file=$outpt
	echo "Outptfile = $outpt"
elif [ "$total_bams" -eq 1 ]; then
	echo "Total files eq 1"
	outpt_file=`echo $bams | sed -e 's/.bam//g'`
elif [ "$total_bams" -gt 1 ]; then
	echo "Total files gt 1"
	errorsms "You have to set outpt"
fi

command="featureCounts ${pair_flag} ${illumina_strand} -t exon -g gene_id -a ${gtf} -o ${outpt_file} -T $procs ${progArgs} $bams"

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export command=$command
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi


#### Example batch run code block

# commandFile="$title.bwa.cmds"
# rm -f $commandFile
# echo -n "" > $commandFile
# while read line
# do
#         l=($line)
#         id=${l[0]}
#         name=${l[1]}
# 
#         if [ $paired -eq 1 ]; then
#             echo "Running in paired end mode : $id"
#             peAlign $id $name > $id.$name.bwa.sh
#         else
#             echo "Running in single end mode"
#             seAlign $id $name > $id.$name.bwa.sh
#         fi
#         #sub_ncores normal bwa 20000 ${procs} $name "sh $id.$name.bwa.sh"
#         echo "${mconf_bashexec} $id.$name.bwa.sh" >> $commandFile
# 
# done < <(grep -v "^#" $inputFile)
# 
# # commands
# 
# farm submission: commandfile var = commandFile
