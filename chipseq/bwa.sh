Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Uses BWA aligner to align reads. 
Gives you sorted bam files.

++ Mandatory args: 
-t|--title [project_Xmen]
-c|--comments [bwa]
-f|--inputFile [input file. see below]
-F|--farm [no input]
-m|--mem [20000]   # becomes bigger with bam file size
-p|--procs [1]
-r|--paired [0/1]
-o|--organism ~/GENOMES/bwa/{genome_name}_bwa


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

mandatory_fails=`mandatory title comments inputFile farm mem procs organism`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi



seAlign() {
    id=$1
    name=$2
    echo "${mconf_bwa} aln -q 15 -t ${procs} $index $id.fq.gz > $name.sai 
    ${mconf_bwa} samse -f $name.sam $index $name.sai $id.fq.gz | ${mconf_samtools} view -@ ${procs} -bS - | ${mconf_samtools} sort -@ ${procs} - -o $name.bam 
    rm -f $name.sai $name.sam"
}


peAlign() {
    id=$1
    name=$2
    echo "${mconf_bwa} aln -q 15 -t ${procs} $index $id.pair1.fq.gz > $name.pair1.sai 
    ${mconf_bwa} aln -q 15 -t ${procs} $index $id.pair2.fq.gz > $name.pair2.sai 
    ${mconf_bwa} sampe $index $name.pair1.sai $name.pair2.sai $id.pair1.fq.gz $id.pair2.fq.gz  | ${mconf_samtools} view -@ ${procs} -bS - | ${mconf_samtools} sort -@ ${procs} - -o ${name}.bam 
    rm -f $name.pair1.sai $name.pair2.sai $name.sam"
}


namefile=$1; shift
index=${organism}
MEM=20000
SCRIPT_PATH=`readlink -f $0`
BIN=`dirname $0`

commandFile="$title.bwa.cmds"
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
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi
