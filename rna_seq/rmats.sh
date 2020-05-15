# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Uses RMATS2 turbo to get differentially spliced events

++ Mandatory args: 
-t|--title [rmats_test]
-f|--inputFile [input file. see below]
-O|--outpt [output folder, overwritten if one exists]
-F|--farm [no input]

++ optional args:
-m|--mem [20000]   # becomes bigger with bam file size
-p|--procs [1]
-P|--progArgs " -t paired " for paired end reads


++ inputFile format
cntrl1.bam	0
cntrl2.bam	0
ko1.bam	1
ko2.bam	1

0 = control samples, 1 = test samples
++ Packages used
rmats2 (turbo version). 

EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title inputFile gtf outpt farm`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

controlfile="48H_scramble_dox.txt"
testfile="120H_dimt1_ko.txt"

tempfile=`${mconf_tempfile_cmd}`

infosms "Creating control file ... ${tempfile}.ctrl.txt"
awk ' $2 == 0  { print $1 } ' ${inputFile} | xargs | sed -e 's/ /,/g' > ${tempfile}.ctrl.txt
awk ' $2 == 1 { print $1 } ' ${inputFile} | xargs | sed -e 's/ /,/g' > ${tempfile}.test.txt


#gtf="/lustre/scratch119/casm/team163gv/SHARED/reference/gtfs/grch38_hg38.gencode.v27.annotation.gtf"
#outpt=`echo ${controlfile}__vs__${testfile} | sed -e 's/.txt//g'`
#otherOptions=" -t paired "

rm -rf ${outpt} ; mkdir ${outpt}
cmd="${mconf_rmatspy} --b1 ${tempfile}.ctrl.txt --b2 ${tempfile}.test.txt --gtf $gtf --od $outpt $progArgs"

echo $cmd



#farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
 	export command=${cmd}
 	$mconf_installdir/bin/farmsub.sh
else 
 	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi


#### Example batch run code block

# commandFile="$title.bwa.cmds"
# rm -f $commandFile
# echo -n "" > $commandFile
# while read line
# do
#         l=($line)
#         id=${l[0]}
#         name=${l[1]}
# 
#         if [ $paired -eq 1 ]; then
#             echo "Running in paired end mode : $id"
#             peAlign $id $name > $id.$name.bwa.sh
#         else
#             echo "Running in single end mode"
#             seAlign $id $name > $id.$name.bwa.sh
#         fi
#         #sub_ncores normal bwa 20000 ${procs} $name "sh $id.$name.bwa.sh"
#         echo "${mconf_bashexec} $id.$name.bwa.sh" >> $commandFile
# 
# done < <(grep -v "^#" $inputFile)
# 
# # commands
# 
# # farm submission: commandfile var = commandFile
# if [ ! -z "$farm" ]; then
# 	export fileOfCommands=$commandFile
# 	$mconf_installdir/bin/farmsub.sh
# else 
# 	warnsms "Not running as farm: Not recommended"
# 	warnsms "$mconf_bashexec $commandFile"
# fi
