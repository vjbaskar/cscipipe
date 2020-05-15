# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-P|--progArgs farm_job_id
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

jid=$progArgs
count=0
sms "Will send email once $jid is done. You can now exit shell"
#count=$((bjobs | awk -v j=$jid ' $1 == j ' | wc -l))
email="$USER@sanger.ac.uk"
#(while `bjobs | awk -v j=$jid ' $1 == j ' | wc -l` ; do sleep 10; done) && mail -s "$jid done; $email"  &

count=1
while [ 1 ]
do
	sleep 500
	count=$((`bjobs | awk -v j=$jid ' $1 == j ' | wc -l`))
	#echo "$jid == $count"
	if [ "$count" == "0" ]; then
		break;
	fi
done
t=`tempfile`
bhist -l $jid > $t
mail -s "HPC job: $jid done" $email < $t



# farm submission: commandfile var = commandFile
# if [ ! -z "$farm" ]; then
# 	export fileOfCommands=$commandFile
# 	$mconf_installdir/bin/farmsub.sh
# else 
# 	warnsms "Not running as farm: Not recommended"
# 	warnsms "$mconf_bashexec $commandFile"
# fi
