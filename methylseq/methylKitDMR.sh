Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Runs methylkit pairwise comparisons

++ Mandatory args: 
-f|--inputFile = input RData file from methylKitCompare 
-0|--cmd0 = distance between methylated bases to merge into DMRs (recommended = 100bp)
-F|--farm = farm option


++ Packages used
R, R:methylKit, R:ChIPSeeker, R:org.*db, R:Txdb*

++ Supported organisms
mm10, hg38

EOM
echo -en "\033[30m"
}

RSCRIPT_NAME="${mconf_installdir}/bin/methylseq/methylKitDMR.R "

input_RData=$inputFile


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory inputFile farm cmd0`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi



base_cmd="${mconf_Rscript} ${RSCRIPT_NAME}"

# if [ ${dataType} == "nondirectional" ]; then
# 	infosms "Turning off directional alignment"
# 	base_cmd="${base_cmd} --non_directional "
# fi


export comments="methylKitDMR"
commandFile="$title.${comments}.cmds"
cmd="${base_cmd} -i $input_RData -w $cmd0"
export command="${cmd}"
# commands
# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$command"
fi
