# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Uses diffBind and a set of peaks to find differential peak binding.

++ Mandatory args: 
-t|--title [project_Xmen]
-f|--inputFile [input file. see below]
-F|--farm [no input]
-m|--mem [20000]   # becomes bigger with bam file size
-p|--procs [1]
-o|--organism [hg19,hg38,mm10]
-1|--cmd1 comparisons file
-O|--outpt output folder


++ inputFile format
SampleID	Condition	Replicate	bamReads	ControlID	bamControl	Peaks	PeakCaller
gRNA_empty1	gRNA_empty	R1	BRD4_Empty_A.bam	IgG_empty1	IgG_Empty_A.bam	BRD4_Empty_A_macs_0.05_peaks_slop250.bed	bed
gRNA_empty2	gRNA_empty	R2	BRD4_Empty_B.bam	IgG_empty2	IgG_Empty_B.bam	BRD4_Empty_B_macs_0.05_peaks_slop250.bed	bed
gRNA_short1	gRNA_short	R1	BRD4_SRPK1_A.bam	IgG_gRNA_short1	IgG_SRPK1_A.bam	BRD4_SRPK1_A_macs_0.05_peaks_slop250.bed	bed
gRNA_short2	gRNA_short	R2	BRD4_SRPK1_B.bam	IgG_gRNA_short2	IgG_SRPK1_B.bam	BRD4_SRPK1_B_macs_0.05_peaks_slop250.bed	bed

++ comparisons file format
gRNA_empty	gRNA_short

++ Packages used
diffBind, R

++ Useful tips
mark duplicates because diffbind requires this and script removes duplicates before counting.


EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title inputFile organism cmd1 outpt`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

#title = Sys.getenv("title")
#sampleFile = Sys.getenv("sampleFile")
#cols = Sys.getenv("comments")
#inputType=Sys.getenv("cmd0")
#outputType=Sys.getenv("cmd1")
#ofold=Sys.getenv("outpt")
comparisons=$cmd1
commandFile="$title.diffBind.cmds"
rm -f $commandFile
echo -n "" > $commandFile
#diffbind="${mconf_Rscript} ${mconf_installdir}/bin/chipseq/diffbind.R"
diffbind="${mconf_R} --no-save --vanilla < ${mconf_installdir}/bin/chipseq/diffbind.R"
while read line
do
	l=($line)
	cond1=${l[0]}
	cond2=${l[1]}
	echo "export title=$title ; export inputFile=$inputFile ; export organism=$organism ; export cond1=$cond1; export cond2=$cond2 ; export outpt=$outpt; $diffbind" >> $commandFile

done < <(grep -v "^#" $comparisons)



# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export comments="diffbind"
	export fileOfCommands=$commandFile
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi



