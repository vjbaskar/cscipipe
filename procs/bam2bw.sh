Help() { 
echo -en "\033[31m"
cat << EOM
DESC:
Convert bam to bigWig

++ Mandatory args: 
	-t|--title
	-c|--comment
	-f|--inputFile
	-o|--organism = hs,mm
	-e|--expType = rnaseq, chipseq
	-F|--farm = flag

++ recommended: 
	-P|--progArgs = optional program specific params [ eg. --normalizeUsing [RPKM, CPM, BPM, RPGC, None], --region CHR:START:END, -r CHR:START:END ]
	Use one of the entered methods to normalize the number of reads per bin. By default, no normalization is performed. RPKM = Reads Per Kilobase per Million mapped reads; CPM = Counts Per Million mapped reads, same as CPM in RNA-seq; BPM = Bins Per Million mapped reads, same as TPM in RNA-seq; RPGC = reads per genomic content (1x normalization); 

	-m|--mem=20000   # becomes bigger with bam file size
	-p|--procs=1

	***** inputFile format
	bamfile1.bam 
	bamfile2.bam 

++ Packages used
	deeptools - https://deeptools.readthedocs.io/en/develop/index.html
	version =  3.2.1


EOM
echo -en "\033[30m"
}

if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile organism expType `

if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

if [ "$genome" == "hs" ]; then
    gsize=3200000000
fi
if [ "$genome" == "mm" ]; then
    gsize=2159570000
fi

binSize=1
minMapQual=1
genome=$organism

commandFile="$title.$comments.cmds"
echo -n "" > $commandFile
while read line
do
	l=($line)
	bam=${l[0]}
	b=`basename $bam .bam`
	cmd=""
	if [ ! -e $bam.bai ]; then
			cmd="samtools index $bam ;"
	fi
	if [ "$expType" == "rnaseq" ]; then
		cmd="$cmd  $mconf_deeptools/bamCoverage -b $bam -o $b.fw.bw  --minMappingQuality $minMapQual -bs $binSize -p $procs --filterRNAstrand forward $progArgs" 
		echo "$cmd" >> $commandFile
		cmd=""
		cmd="$cmd  $mconf_deeptools/bamCoverage -b $bam -o $b.rev.bw  --minMappingQuality $minMapQual -bs $binSize --normalizeUsing RPKM -p $procs --filterRNAstrand reverse  $progArgs" 
        echo "$cmd" >> $commandFile
    else
        cmd="$cmd $mconf_deeptools/bamCoverage -b $bam -o $b.bw  --minMappingQuality $minMapQual -bs $binSize --normalizeUsing RPKM -p $procs  $progArgs" 
        echo "$cmd" >> $commandFile
    fi
done < $inputFile

#farm submission
if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	export concurrentJobs=20
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi
