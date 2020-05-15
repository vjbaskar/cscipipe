# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM

DESC:
Farm history in the present dir 

++ Mandatory args:
-P|--progArgs [pretty,csv]

++ Optional args:
-c|--comments [search string]

EOM
echo -en "\033[30m"
}

if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory progArgs`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

$mconf_conda_py3
if [ ! -z "$comments" ]; then
	comments="--comments $comments"
fi
python $mconf_installdir/bin/farm/farmhist.py --progArgs $progArgs $comments


# farm submission: commandfile var = commandFile
# if [ ! -z "$farm" ]; then
# 	export fileOfCommands=$commandFile
# 	$mconf_installdir/bin/farmsub.sh
# else 
# 	warnsms "Not running as farm: Not recommended"
# 	warnsms "$mconf_bashexec $commandFile"
# fi


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
