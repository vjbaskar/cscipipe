# Get data from irods

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-s|--sampleFile = sample file with specific format
-p|--procs = total parallel downloads
-0|--cmd0 = file type [tsv,csv]

**** Packages used
irods

EOM
echo -en "\033[30m"
}

if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title sampleFile procs cmd0`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

#Set Script Name variable
THREADS=$procs



function md5sumCheck(){
	irods_file=$1
	local_file=$2
	returnVal=0
	md5=`imeta ls -d $irods_file md5 | awk ' $1 ~ /value/ { print $2 } '`
	lmd5=`md5sum $local_file ! awk ' { print $1 } '`
	echo $md5
	echo $lmd5
	if [ $md5 == $lmd5 ]; then
		returnVal=1
	fi
	echo $returnVal
}

# Get data from irods
function getDataFromIrods(){

function md5sumCheck(){
	irods_file=$1
	local_file=$2
	returnVal=0
	md5=`imeta ls -d $irods_file md5 | awk ' $1 ~ /value/ { print $2 } '`
	lmd5=`md5sum $local_file | awk ' { print $1 } '`
	if [ $md5 == $lmd5 ]; then
		returnVal=1
	fi
	echo $returnVal
}

	filetype="cram"
	runId=$1; shift
	sample=$1; shift # MOLM_150_C
	lanes=$1; shift # 3,4,5,6
	lanes=`echo "$lanes" | sed -e "s/,/ /g"`

	sms " runid = $runId | sample name = $sample | file type = $filetype | lanes = $lanes" 

	ofile=${sample}
	targetFolder=$sample
	
	rm -rf $targetFolder
	if [ ! -z $targetFolder ]; then
		sms "Dumping data in $targetFolder "
		mkdir $targetFolder
	fi
#imeta qu -z seq -d id_run = $id_run and target = 1 and manual_qc = 1 and sample = $sample > $imetafile
	rm -f $imetafile

	mkdir ${targetFolder}/finalfiles/
	for i in $lanes
		do
			md5checkVal=0
			imetafile=$targetFolder/${runId}.${sample}.$i.imeta
			echo -en "\r \033[33m>> Searching lane, $i \033[0m"
			imeta_cmd="imeta qu -z seq -d type = $filetype and target = 1 and lane = $i and id_run = $runId and sample = $sample"
			$imeta_cmd > $imetafile
			echo $imeta_cmd >> $cwd/$title.irods.cmds
			wc=`cat $imetafile | grep cram | wc -l`
			p=`pwd`
			if [ $wc -gt 0 ]; then
				seqfolder=(`awk ' $1 ~ "collection:" { print $2 } ' $imetafile`)
				files=(`awk ' $1 ~ "dataObj:" { print $2 } ' $imetafile`)
				fileName=${seqfolder}/${files}
				mkdir -p $targetFolder/${runId}
				iget $fileName $targetFolder/${runId}/${files}
				md5checkVal=$(md5sumCheck $fileName $targetFolder/${runId}/${files})
				if [ "$md5checkVal" -eq 1 ]; then
					mv $targetFolder/${runId}/${files} $targetFolder/finalfiles/
				else
					echo "md5checksum failed"
				fi
				echo "$runId	$sample	$i	$seqfolder/$files	mdcheck=$md5checkVal" >> $title.irods.log
					
			fi
		done
		echo ""
	return 0
}


#### codes start
#### Get data
#$mconf_Rscript $mconf_installdir/bin/misc/createProject.R --sa"pleSheet $sampleFile -c RUNID,SAMPLEID,LANES --tsv 1 | sed -e '1d' >| samples.irods.txt

export cmd1="tsv" 
#export outpt=`echo $sampleFile | sed -e "s/$cmd0/irods.txt/" `
export outpt="$title.irods.txt"
export comments="RUNID,SAMPLEID,LANES"
#echo "$title $sampleFile $inputType $cmd0 $cmd1 $outpt $comments" 


$mconf_Rscript $mconf_installdir/bin/procs/filterCols.R 
rm -f $title.irods.log
export -f getDataFromIrods
#getDataFromIrods 24076 3960STDY7100225     1
rm -f $title.irods.log
cat $title.irods.txt | grep -v RUNID | xargs -n3 -I@ -P $THREADS bash -c "$(declare -f getDataFromIrods); getDataFromIrods @"

