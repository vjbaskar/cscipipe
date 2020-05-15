# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Converts bam to bedgraph

++ Mandatory args: 
-t|--title
-c|--comment
-f|--inputFile
-g|--genome = hs,mm
-e|--expType = rnaseq, chipseq
-F|--farm = flag
-m|--mem=20000   # becomes bigger with bam file size
-p|--procs=1


***** inputFile format
bamfile1.bam 
bamfile2.bam 


++ Packages used
samtools, deepTools

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

if [ "$genome" == "hs" ]; then
    gsize=3200000000
fi
if [ "$genome" == "mm" ]; then
    gsize=2159570000
fi

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
	commonCommand="samtools index $bam ; $mconf_deeptools/bamCoverage --minMappingQuality 1 -bs 2 --normalizeUsingRPKM -p $procs -of bedgraph -v -b $bam -o $b.fw.bdg "
	if [ "$expType" == "rnaseq" ]; then
		#cmd="$cmd  $mconf_deeptools/bamCoverage -b $bam -o $b.fw.bdg  --minMappingQuality 1 -bs 2 --normalizeUsingRPKM -p $procs --filterRNAstrand forward -of bedgraph -v True" 
		cmd="$commonCommand --filterRNAstrand forward "
		echo "$cmd" >> $commandFile
		cmd=""
		#cmd="$cmd  $mconf_deeptools/bamCoverage -b $bam -o $b.rev.bdg  --minMappingQuality 1 -bs 2 --normalizeUsingRPKM -p $procs --filterRNAstrand reverse -of bedgraph -v True" 
		cmd="$commonCommand --filterRNAstrand reverse "
        echo "$cmd" >> $commandFile
    else
        #cmd="$cmd $mconf_deeptools/bamCoverage -b $bam -o $b.bdg  --minMappingQuality 1 -bs 2 --normalizeUsingRPKM -p $procs -of bedgraph -v True" 
        cmd=$commonCommand
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
