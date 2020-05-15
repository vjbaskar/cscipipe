# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments
-f|--inputFile
-s|--sampleFile = miso config file

++ recommended: 
-F|--farm = flag
-m|--mem=20000   # becomes bigger with bam file size
-p|--procs=1

	===== Input config file =====
	indexed_events=~/GENOMES/MISO/GRCm38/all_events
	genome_gff3=~/GENOMES/SANGER/GRCm38/Mus_musculus.GRCm38.90.gff3
	control_bam=NPM1_WT.merge.bam
	control_name=NPM1_WT
	test_bam=NPM1_MET16.merge.bam
	test_name=NPM1_MET16
	read_length=75
	frag_length=170
	frag_length_stdev=48
	output_dir=NPM1_miso
	miso_settings_file=miso.settings


	===== Miso config file =====
	[data]
	filter_results = True
	min_event_reads = 20
	strand = fr-firststrand
	[cluster]
	cluster_command = bsub -M 10000 -R "select[mem>10000] rusage[mem=10000]" -G team163-grp
	long_queue_name = normal
	short_queue_name = small
	[sampler]
	burn_in = 500
	lag = 10
	num_iters = 5000
	num_processors = 1

	===== steps =====
	
	fragLength = computes fragment length in RNA-seq data for control and test
	miso = runs the main miso program. Does not require to be submitted by sub..farm because internally it can submit over cluster using settings from miso.settings
	log = optional collection of logs. miso submits chunks of genes as jobs and this command collects all the log files. **** Use to check if miso run is complete ****
	summarise = summarises miso runs for each condition
	compare = compares the summaries of the two conditions
	filter = filters based on certain cut offs
	annotate = annotates the filtered events based on data from indexed_events

	===== Output data =====
	Step	|	Output
	------------------
	miso	|	outputdir/condn/scripts_output # eg. WT_MET16_miso/WT/scripts_output

**** Packages used
miso
https://miso.readthedocs.io/en/fastmiso/

EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

#### Code starts here

getMeanStd(){
    bamfile=$1
    c=`head -1 ${output_dir}/insert-dist/${bamfile}.insert_len`
    declare -a params
    params=(`echo $c | tr "," "\n"`)
    mean=`echo ${params[0]} | awk -F "=" ' { print $2 }'`
    stdev=`echo ${params[1]} | awk -F "=" ' { print $2 }'`
    echo "$mean $stdev"
}


#### Input data
# cfg file
#cfgFile=$inputFile
#sampleFile=$sampleFile
com=${comments}
#output_dir=${title}_${comments}


# Step you want to run
step=$1; shift
R_SCRIPTS="$mconf_installdir/rmds/"
p=`pwd`

# load config
source $p/$inputFile

# Create output directory
mkdir -p ${output_dir}/

# Name you want to use
runName="${control_name}.${test_name}"

sms "Run Name = $runName"
# Get fragment lengths

step="fraglength"
echo "STEP: $step"
t=`basename ${genome_gff3} .gff3`
exons=${t}.min_1000.const_exons.gff
echo "exon_utils --get-const-exons ${genome_gff3} --min-exon-size 1000 --output-dir ${output_dir}/exons/" > ${output_dir}/fraglength.sh
for bamfile in ${control_bam} ${test_bam}
do
	echo "pe_utils --compute-insert-len $bamfile ${output_dir}/exons/$exons --output-dir ${output_dir}/insert-dist/"
done >> ${output_dir}/fraglength.sh


export comments="${com}.${step}"
export queue="normal"
export mem="10000"
export command="sh ${output_dir}/fraglength.sh"  

$mconf_bashexec $mconf_installdir/bin/farmsub.sh

export jobIDtoWait="${title}.${comments}"

#echo $jobIDtoWait


# Run MISO
#if [ "$step" == "miso" ]; then
step="miso"
echo "STEP: $step"
#controlPE=`getMeanStd ${control_bam}`
#testPE=`getMeanStd ${test_bam}`
# ${frag_length} ${frag_length_stdev}
echo "
getMeanStd(){
    bamfile=\$1
    c=\`head -1 ${output_dir}/insert-dist/${bamfile}.insert_len\`
    declare -a params
    params=(\`echo \$c | tr \",\" \"\n\"\`)
    mean=\`echo \${params[0]} | awk -F "=" ' { print \$2 }'\`
    stdev=\`echo \${params[1]} | awk -F "=" ' { print \$2 }'\`
    echo \"\$mean \$stdev\"
}
controlPE=\`getMeanStd ${control_bam}\`
testPE=\`getMeanStd ${test_bam}\`
miso --run ${indexed_events}/  ${control_bam} --output-dir ${output_dir}/${control_name}  --read-len ${read_length} --paired-end ${controlPE} --settings-file ${miso_settings_file} --use-cluster --no-wait --chunk-jobs 500
miso --run ${indexed_events}/  ${test_bam} --output-dir ${output_dir}/${test_name}  --read-len ${read_length} --paired-end ${testPE} --settings-file ${miso_settings_file} --use-cluster --no-wait --chunk-jobs 500
" > ${output_dir}/miso.sh
export comments="${com}.${step}"
export queue="normal"
export mem="500"
export command="sh ${output_dir}/miso.sh"  
$mconf_installdir/bin/farmsub.sh
unset jobIDtoWait


