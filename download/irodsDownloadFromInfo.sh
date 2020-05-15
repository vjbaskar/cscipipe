# Get data from irods


if [ $args_return  -gt 0 ]; then
print_stderr "
++ Mandatory args: 
-t|--title
-s|--sampleFile = sample file. A file generated from irodsRunInfo or irodsSampleInfo
-p|--procs = total parallel downloads
"
exit 2
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

	fileName=$2
	runId=$3
	sample=$1
	targetFolder=${sample}
	files=`basename $fileName`

	mkdir -p ${targetFolder}/finalfiles/
	mkdir -p ${targetFolder}/${runId}/
	
	echo "file = $fileName"
	echo "id = $runId"
	echo "sample = $sample"
	echo "to = $files"
	echo "-"
	
 	md5checkVal=0
 	iget $fileName $targetFolder/${runId}/${files}
 	md5checkVal=$(md5sumCheck $fileName $targetFolder/${runId}/${files})
	if [ "$md5checkVal" -eq 1 ]; then
		mv $targetFolder/${runId}/${files} $targetFolder/finalfiles/
	else
		echo "md5checksum failed"
	fi
	echo "$runId	$sample	$i	$seqfolder/$files	mdcheck=$md5checkVal" >> ${title}.irods.log
	echo ""
	return 0
}


#### codes start
#### Get data

cat $sampleFile | awk ' { print $1, $2, $3 } ' | grep -v "sample_id" | xargs -n3 -I@ -P $THREADS bash -c "$(declare -f getDataFromIrods); getDataFromIrods @"
