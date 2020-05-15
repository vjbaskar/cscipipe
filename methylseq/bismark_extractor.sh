Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Uses BWA aligner to align reads. 
Gives you sorted bam files.

++ Mandatory args: 
-t|--title [project_Xmen]
-c|--comments [bismark_methyl]
-f|--inputFile [input file. see below]
-F|--farm [no input]
-m|--mem [10000]   # becomes bigger with bam file size
-p|--procs [ 1 ]
-r|--paired [0/1]
-P|--progArgs [optional for providing with other flags right out of the manual]
-o|--organism [ repository/bismark/indexfolder. Should contain genome.fa ]

++ inputFile format (eg. A1_DMSO1.pair1.fq.gz and A1_DMSO.pair2.fq.gz
DP_E2_OX_bismark_outpt/DP_E2_OX.pair1_bismark_bt2_pe.bam
DP_E4_BS_bismark_outpt/DP_E4_BS.pair1_bismark_bt2_pe.bam
DP_F5_OX_bismark_outpt/DP_F5_OX.pair1_bismark_bt2_pe.bam
DP_F6_BS_bismark_outpt/DP_F6_BS.pair1_bismark_bt2_pe.bam
DP_G12_BS_bismark_outpt/DP_G12_BS.pair1_bismark_bt2_pe.bam
DP_G7_OX_bismark_outpt/DP_G7_OX.pair1_bismark_bt2_pe.bam
DP_I13_OX_bismark_outpt/DP_I13_OX.pair1_bismark_bt2_pe.bam
DP_I14_BS_bismark_outpt/DP_I14_BS.pair1_bismark_bt2_pe.bam

++ Packages used
samtools, bismark_methylation_extractor


EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile farm mem procs paired organism`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi



seAlign() {
    id=$1
    ofold=$2
    echo "${base_cmd} -s -o ${ofold} $id "
}


peAlign() {
    id=$1
    ofold=$2
    echo "${base_cmd} -p -o ${ofold} $id "
}

# However, please note that a typical Bismark run will use several cores already (Bismark itself, 2 or 4 threads of
# Bowtie/Bowtie2, Samtools, gzip etc...) and ~10-16GB of memory depending on the choice of aligner
# and genome. WARNING: Bismark Parallel (BP?) is resource hungry! Each value of --parallel specified
# will effectively lead to a linear increase in compute and memory requirements, so --parallel 4 for
# e.g. the GRCm38 mouse genome will probably use ~20 cores and eat ~40GB or RAM, but at the same time
# reduce the alignment time to ~25-30%. You have been warned.

nprocs=`echo $procs/8 | bc`
#procs=1
base_cmd="${mconf_bismark}/bismark_methylation_extractor --bedGraph --parallel $procs --cytosine_report --genome_folder $organism $progArgs "

# if [ ${dataType} == "nondirectional" ]; then
# 	infosms "Turning off directional alignment"
# 	base_cmd="${base_cmd} --non_directional "
# fi



commandFile="$title.${comments}.cmds"
rm -f $commandFile
echo -n "" > $commandFile
while read line
do
        l=($line)
        id=${l[0]}
        ofold=`echo $id | sed -e 's/bam$/meOutpt/g'`

        if [ $paired -eq 1 ]; then
            echo "Running in paired end mode : $id"
            cmd=`peAlign $id $ofold`
        else
            echo "Running in single end mode"
            cmd=`seAlign $id $ofold`
        fi
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
