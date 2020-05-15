
Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments

++ recommended: 
-q|--queue=normal
-m|--mem=20000   # becomes bigger with bam file size
-p|--procs=1
-x|--command

++ Look above in ARGUMENTS for more options

EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi



## FUnctions for setting the variables if they are optional
optional(){
	if [  -z "$queue" ]; then
		warnsms "Queue is undefined. Setting it to normal"
		export queue="normal"
	fi
	if [  -z "$mem" ]; then
		warnsms "Memory not set. Setting it to default = 1000 mb"
		export mem=1000
	fi
	if [  -z "$procs" ]; then
		warnsms "procs not set. Setting it to default = 1"
		export procs=1
	fi

	if [ ! -z "${fileOfCommands}" ]; then
		if [  -z "$concurrent" ]; then
			warnsms "concurrent jobs not set. Setting it to default = 20"
			export concurrent=20
		fi
	fi
}

mandatory(){
	echo "Mandatory Args"
}

## Functions to generate bsub script

# farmGroup(){
# 	c=`echo $HOSTNAME | grep cgp -c`
# 	if [ "$c" -eq "0" ]; then
# 		echo "#BSUB -G ${farmgroup}"
# 	fi
# }

coreBSUB() {
# - core bsub function
	#command=$1
	d=`date +%_F.%Hhr`
	ofile="${title}.${comments}"
echo "
#!${mconf_bashexec}
#source ~/.bashrc
#BSUB -M ${mem} 
#BSUB -R \"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1] \" 
#BSUB -n $procs 
#BSUB -q ${queue}"
##BSUB -B"
#farmGroup
}

jobName(){
# - Job name
	echo "#BSUB -J ${ofile}"
}

bsubFiles(){
# - Out and err files
	echo "
#BSUB -o .bsub/${ofile}.farm
#BSUB -e .bsub/${ofile}.farm"
}
waitBSUB(){
# - Wait till finish
	echo "#BSUB -K"
}

waitForJob(){
# - Wait till finish of a specific job
echo "#BSUB -w 'done(${jobIDtoWait})'"
}

arrayBSUB(){
# - Run as array
	totalJobs=`cat $fileOfCommands | wc -l`
	if [ -z ${concurrent} ]; then
		concurrent=${totalJobs}
	fi
echo "
#BSUB -o .bsub/${ofile}.%I.farm
#BSUB -e .bsub/${ofile}.%I.farm
#BSUB -J \"${ofile}[1-$totalJobs]%${concurrent}\"
command=\`sed -e \"\${LSB_JOBINDEX}q;d\" $fileOfCommands\`
echo \"\$command\"
eval \"\$command\"

"
#echo \"\$command\" | sed -e 's/;/\n/g' | ${mconf_parallel}/parallel -j ${procs} eval
}

timeBSUB(){
	t=$1
echo "#BSUB -W $t"
}


commandBSUB(){
echo "
set -euxo pipefail
echo \"Job id = \$LSB_JOBID \"
d=\`date\`
echo \"[ \$d ] ${ofile} starts\"
$command
d=\`date\`
echo \"[ \$d ] ${ofile} ends\"
" 
}

writeBSUB() {
	# -- A function to create submission scripts for LSF

	if [ ! -d .bsub ]; then
			mkdir .bsub
	fi

	coreBSUB
	if [ ! -z "$wait_to_finish" ]; then
		waitBSUB
	fi
	if [ ! -z "$jobIDtoWait" ]; then
		warnsms "Requesting to wait for $jobIDtoWait to finish"		
		waitForJob
	fi
	if [ ! -z "$walltime" ]; then
		timeBSUB $walltime
	fi
	if [ ! -z "$fileOfCommands" ]; then
		command=`arrayBSUB`
	else
		bsubFiles
		jobName
	fi
commandBSUB

}



#### End of functions

# If the values are undefined, default them
optional

# Generate the bsub script
writeBSUB "$command" > ${title}.${comments}.bsub

# Submit job
temp=`$mconf_tempfile_cmd`

#echo "bsub < ${title}.${comments}.bsub > $temp"
bsub < ${title}.${comments}.bsub > $temp
# Get job id
jid=`head -1 $temp | awk ' { print $2 } ' | sed -e "s/>//g" | sed -e "s/<//g"`
infosms ">> Job title: ${title}.${comments}.bsub >> Job ID = $jid"
echo "[date = `date`] [operation = $title] [comments = $comments] [bsub_file = $ofile] [job_id = $jid] [log = ${ofile}] [command = ${command}]" >> joblists.lsf

#echo "" >> joblists.lsf
echo "$jid"
rm -f $temp

