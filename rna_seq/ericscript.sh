# Script: EricScript for fusion detection
# Note that this version needs samtools < 0.19
# Install modules to switch to samtools 0.19 or change it in your bashrc to default to samtools version < 0.19

if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: 
	t|title
	c|comments
	f|inputFile
	0|file0 = database file in $mconf_ericscript_db in cscipipe.conf
	
	++ recommended: 
	F|farm
	m|em=10000   
	p|procs=1
	q|queue=normal

	***** inputFile format
	fastqBaseNameExpt1
	fastqBaseNameExpt2

"
	exit 2
fi

p=`pwd`
loadSamtools="module remove hgi/samtools/0.1.19; module load hgi/samtools/0.1.19"

cp -R $mconf_ericscript $p
ericscript=`basename $mconf_ericscript`

commandFile="$title.$comments.cmds"
echo -n "" > $commandFile
while read line
do
	l=($line)
	id=${l[0]}
	cmd="./$ericscript/ericscript.pl -db $file0 -name $id -o ${id}_ericscript ${id}.*.fq.gz -p $procs"
	echo $cmd
	cmd="$loadSamtools ; $cmd"
	echo $cmd >> $commandFile
done < $inputFile	
#sub_ncores normal star $mem $threads $name "sh $id.$name.star.sh"

if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	export concurrentJobs=20
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi
