Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Runs methylkit pairwise comparisons

++ Mandatory args: 
-t|--title [project_Xmen]
-f|--inputFile [input file. see below]
-F|--farm [no input]
-m|--mem [10000]   # becomes bigger with bam file size
-o|--organism [ mm10,hg38 ]
-O|--outpt [ output prefix ]
-v|--qval [ qvalue cutoff. 0.05 {0.0-1.0} ]
-0|--cmd0 [ percentage diff methylation cutoff. 10 { 0-100 } ]
-1|--cmd1 [ low count cutoff. 10 {0-inf} ]
-2|--cmd2 [ Tiling Window Size. 200 {0-inf} ]
-3|--cmd3 [ Tiling Step Size. 50 {0-inf} ]
-4|--cmd4 [ CpG Bed file. Download from UCSC Table Browser. 4 Fields ]
-p|--procs [ total processors ]

++ inputFile format (bam files should be sorted and indexed)
HFMX7_1E_OX.pair1_bismark_bt2_pe_pba/HFMX7_1E_OX.pair1_bismark_bt2_pe_CpG.txt	HFMX7_1E_OX	0
HFMX7_1F_OX.pair1_bismark_bt2_pe_pba/HFMX7_1F_OX.pair1_bismark_bt2_pe_CpG.txt	HFMX7_1F_OX	1
HFMX7_1H_OX.pair1_bismark_bt2_pe_pba/HFMX7_1H_OX.pair1_bismark_bt2_pe_CpG.txt	HFMX7_1H_OX	1
HFMX7_1I_OX.pair1_bismark_bt2_pe_pba/HFMX7_1I_OX.pair1_bismark_bt2_pe_CpG.txt	HFMX7_1I_OX	0

Notes: 0 = control; 1 = test

++ Packages used
R, methylKit, genomation


EOM
echo -en "\033[30m"
}

RSCRIPT_NAME="${mconf_installdir}/bin/methylseq/methylKitCompare.R "

outpt_prefix=$outpt
sampleFile=$inputFile
qcut=$qval
percCutOff=$cmd0
lowcountCutoff=$cmd1
tilingWindowSize=$cmd2
tilingStepSize=$cmd3
cpgBedFile=$cmd4


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title inputFile`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi



base_cmd="${mconf_Rscript} ${RSCRIPT_NAME}"

# if [ ${dataType} == "nondirectional" ]; then
# 	infosms "Turning off directional alignment"
# 	base_cmd="${base_cmd} --non_directional "
# fi


export comments="methylcomp"
commandFile="$title.${comments}.cmds"
cmd="${base_cmd} -O $outpt_prefix -f $sampleFile -v $qcut -c $percCutOff -o $organism -p $procs -l $lowcountCutoff -w $tilingWindowSize -s $tilingStepSize -b $cpgBedFile"
export command="${cmd}"
# commands
# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$command"
fi
