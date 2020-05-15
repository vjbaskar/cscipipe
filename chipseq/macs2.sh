Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title 
-f|--inputFile : inputfile
-r|--paired [0,1]
-v|--qval [0-1]
-m|--mem [10000]
-p|--procs [1]
-o|--organism [hs,mm]
-c|--comments [broadPeak,narrowPeak]
-F|--farm 

***** inputFile format
bamfile1.bam control1.bam	outpt1
bamfile2.bam control1.bam	outpt2


**** Packages used
macs2 

EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title organism paired qval inputFile mem procs comments farm`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

# Conda env installation of macs2
#source activate p27

commandFile="$title.macs.cmds"
echo $commandFile
rm -f $commandFile


cmd_callpeak="${mconf_macs2_exec} callpeak  "
cmd_createbdg="-B --SPMR "
cmd_broadpeak="--broad --broad-cutoff 0.1 "
cmd_verybroadpeak="--broad --broad-cutoff 0.25 "
while read line
	do
	l=($line)
	bamfile=${l[0]}
	controlfile=${l[1]}
	outptPrefix=${l[2]}
	sms $outptPrefix
	if [ $paired -eq 1 ]; then
		fileFormat="BAMPE"
	else
		fileFormat="BAM"
	fi
	qval=0.05
	binit=`basename $bamfile .bam` 
	b=`echo "${b}_macs_${qval}"`
	c=`basename $controlfile .bam`
	comments=`toLower $comments`
#	odir="bw"
#	macs_bw="${mconf_macs2}/macs2 callpeak -t $bamfile -c $controlfile -f $fileFormat -g $organism -n $b -q 0.99 -B --SPMR --outdir $odir"
#	echo $macs_np >> $commandFile
	cmd_comparison="-t $bamfile -c $controlfile -f $fileFormat -g $organism -n $binit "
	macs_bdg="$cmd_callpeak $cmd_comparison $cmd_createbdg -q $qval --outdir bw"
	echo $macs_bdg >> $commandFile
	for qval in 0.2 0.1 0.05 0.01 0.00001
	do
		b=`echo "${binit}_macs_${qval}"`
		cmd_comparison="-t $bamfile -c $controlfile -f $fileFormat -g $organism -n $b "
		if [ "$comments" == "narrowpeak" ]; then
			#sms "Running in narrow peak mode"
			odir="${outptPrefix}.np"
			macs_np="$cmd_callpeak $cmd_comparison -q $qval --outdir $odir"
			echo $macs_np >> $commandFile
		fi
		if [ "$comments" == "broadpeak" ]; then
			sms "Running in broad peak mode"		
			odir="${outptPrefix}.bp"
			macs_bp="$cmd_callpeak $cmd_comparison -q $qval --outdir $odir $cmd_broadpeak"
			odir="${outptPrefix}.vbp"
			macs_vbp="$cmd_callpeak $cmd_comparison -q $qval --outdir $odir $cmd_verybroadpeak"
			echo $macs_bp >> $commandFile
		fi
	done
	
done < $inputFile

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export comments="macs2"
	export fileOfCommands=$commandFile
 	$mconf_installdir/bin/farmsub.sh
else 
 	warnsms "Not running as farm: Not recommended"
 	warnsms "$mconf_bashexec $commandFile"
fi
