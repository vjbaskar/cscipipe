Help() { 
echo -en "\033[31m"
cat << EOM
	BigWig file for chipseq from bedgraph file
	++ Mandatory args: 
	-t|--title 
	-f|--inputFile : inputfile
	-y|--chrType : chromosome type [ucsc, ensembl]
	-m|--mem [10000]
	-p|--procs [1]
	-o|--organism [/lustre/scratch119/realdata/mdt1/team163/SHARED/reference/genomeSizes/{genomeName}.chr.sizes]
	-F|--farm 

	==== inputfile
	10_MUT_NPM1_macs_0.05.spmr_treat_pileup.bdg
	10_TY1_macs_0.05.spmr_treat_pileup.bdg
	10_WT_NPM1_macs_0.05.spmr_treat_pileup.bdg
	30_MUT_NPM1_macs_0.05.spmr_treat_pileup.bdg
	30_TY1_macs_0.05.spmr_treat_pileup.bdg
	30_WT_NPM1_macs_0.05.spmr_treat_pileup.bdg

	****
	UCSC tools: bedGraphToBigWig
	Bedtools: sortBed

EOM
echo -en "\033[30m"2
}

if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title inputFile chrType mem procs organism farm`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

#genome=$organism
gtype=$chrType
bedGraphList=$inputFile
genomeSizes=$organism
#mem=$mem
#nProcs=$procs

commandFile="$title.bdg2bw.cmds"
rm -f $commandFile

while read line
do
	l=($line)
	file=${l[0]}
	echo "File = $file"
	f=`basename $file .bdg`
	#tfile=`tempfile -d ./`
	tfile=`${mconf_tempfile_cmd}`
	command=""
	
	if [ "$gtype" == "ensembl" ]; then
		command=" awk ' OFS=\"\t\" { if(\$1 ~/^[0-9]|^X|^Y/) print \"chr\"\$0 }  ' $file | sortBed -i - > $tfile"
	else 
		command="sortBed -i $file > $tfile"
	fi
	 command="$command ; bedGraphToBigWig $tfile $genomeSizes $f.bw ; rm -f $tfile" 
	 echo $command >> $commandFile

done < $bedGraphList

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export comments="bdg2bw"
	export fileOfCommands=$commandFile
	echo $commandFile
 	$mconf_installdir/bin/farmsub.sh
else 
 	warnsms "Not running as farm: Not recommended"
 	warnsms "$mconf_bashexec $commandFile"
fi
