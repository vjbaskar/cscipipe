Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Uses BWA aligner to align reads. 
Gives you sorted bam files.

++ Mandatory args: 
-t|--title [project_Xmen]
-c|--comments [bismark]
-f|--inputFile [input file. see below]
-d|--dataType [directional, nondirectional]
-F|--farm [no input]
-m|--mem [50000]   # becomes bigger with bam file size
-p|--procs [ > 8, in multiples of 8 because approximately bismark 1-core ~= farm 8-core ]
-r|--paired [0/1]
-o|--organism /lustre/scratch119/casm/team163gv/SHARED/reference/bismark/* [mm10]
-P|--progArgs [optional for providing with other flags: eg. --bowtie1]


++ inputFile format (eg. A1_DMSO1.pair1.fq.gz and A1_DMSO.pair2.fq.gz
A1_DMSO1
A2_DMSO2
A3_DRUG1
A4_DRUG2

++ Packages used
samtools, bismark, bowtie1 or 2


EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile mem procs organism dataType`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi



seAlign() {
    id=$1
    lcmd=$2
    echo "${lcmd} ${id}.fastq.gz"
}


peAlign() {
    id=$1
    lcmd=$2
    echo "${lcmd} -1 ${id}.pair1.fq.gz -2 ${id}.pair2.fq.gz "
}

# However, please note that a typical Bismark run will use several cores already (Bismark itself, 2 or 4 threads of
# Bowtie/Bowtie2, Samtools, gzip etc...) and ~10-16GB of memory depending on the choice of aligner
# and genome. WARNING: Bismark Parallel (BP?) is resource hungry! Each value of --parallel specified
# will effectively lead to a linear increase in compute and memory requirements, so --parallel 4 for
# e.g. the GRCm38 mouse genome will probably use ~20 cores and eat ~40GB or RAM, but at the same time
# reduce the alignment time to ~25-30%. You have been warned.

nprocs=`echo $procs/8 | bc`
#procs=1
base_cmd="${mconf_bismark}/bismark --genome $organism --parallel $nprocs "

if [ ${dataType} == "nondirectional" ]; then
	infosms "Turning off directional alignment"
	base_cmd="${base_cmd} --non_directional "
fi



commandFile="$title.${comments}.cmds"
rm -f $commandFile
echo -n "" > $commandFile
while read line
do
        l=($line)
        id=${l[0]}
        cmd="${base_cmd} -o ${id}_bismark_outpt --temp_dir ${id}_bismark_tempdir "
	if [ $paired -eq 1 ]; then
            echo "Running in paired end mode : $id"
            cmd=`peAlign $id "$cmd"`
        else
            echo "Running in single end mode"
            cmd=`seAlign $id "$cmd"`
        fi
        echo "${cmd}" >> $commandFile
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
