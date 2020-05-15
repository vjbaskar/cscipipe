# All args come from cscipipe command

export nohistory=1 ; 
Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-P|--progArgs: cat, more, less, edit, grep


++ optional: 
-c|--comments: simple search. For eg. use date, command names etc. Uses linux grep.
-i|--id: local

++ example:

cscipipe hist -P cat
cscipipe hist -P cat -c "Wed Apr 2"

EOM
echo -en "\033[30m"
}

if [ "$id" == "local" ]; then
	cscipipe_cmds="./.cscipipe/hist.cmds"
else
	cscipipe_cmds="$HOME/.cscipipe/hist.cmds"
fi



if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory progArgs`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

if [ ! -z "$comments" ]; then
	grep "$comments" $cscipipe_cmds
else 
	if [ "$progArgs" == "edit" ]; then
		progArgs="vim"
	fi
	if [ "$progArgs" == "grep" ]; then
		progArgs="grep -i "
	fi
	$progArgs $cscipipe_cmds
fi
