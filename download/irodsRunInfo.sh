# Gets sample information given a set of sample ids
# All args come from cscipipe command

if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: 
	-t|--title
	-f|--inputFile	

	***** inputFile format
	24847
	21990
"
	exit 2
fi

function getType(){
	
	tempfile=$1
	dataType=$2
	w=`cat $tempfile | grep -A 1 "$dataType$" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
	echo $w
}

function getMetaData(){
	
#getMetaData /seq/21908/21908_1#1.cram
    fileName=$1
    t=`tempfile -d ./`
    imeta ls -d $fileName > $t
    
   # fn=`cat $t | grep -A 1 "sample_supplier_name" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    fn=`getType $t "sample_supplier_name"`
    lane=`getType $t "lane"`
    tagIndex=`getType $t "tag_index"`
    sample_supplier_name=`getType $t "sample_supplier_name"`
    seqchksum=`getType $t "seqchksum"`
    total_reads=`getType $t "total_reads"`
    library_type=`getType $t "library_type"`
    sample_donor_id=`getType $t "sample_donor_id"`
    reference=`getType $t "reference"`
    id_run=`getType $t "id_run"`
    is_paired_read=`getType $t "is_paired_read"`
    printf "sample_id\tfile_name\tid_run\tsample_donor_id\tsample_supplier_name\treference\tlibrary_type\tis_paired_read\n"
    printf "${sample}\t$fileName\t$id_run\t$sample_donor_id\t$sample_supplier_name\t$reference\t$library_type\t$is_paired_read\n"
    rm -f $t
    
}

filetype="cram"
counts=0
while read line
do
	l=($line)
	runid=${l[0]}
	sms "Retrieving data from irods >> $runid " 

	ofile=${runid}
	
	targetFolder=$runid
	rm -rf targetFolder
	imetafile=${ofile}.imeta
	#imeta qu -z seq -d id_run = $id_run and target = 1 and manual_qc = 1 and sample = $sample > $imetafile
	#imeta qu -z seq -d  target = 1 and manual_qc = 1 and sample = $sample > $imetafile
	imeta qu -z seq -d id_run = $runid > $imetafile
	
	wc=`cat $imetafile | wc -l`

	#getMetaData /seq/24684/24684_2#1.cram

	#### If the imetafile contains some lines then do the following
	p=`pwd`

	if [ $wc -gt 0 ]; then
		seqfolder=(`awk ' $1 ~ "collection:" { print $2 } ' $imetafile`)
		files=(`awk ' $1 ~ "dataObj:" { print $2 } ' $imetafile`)
		totalfiles=`expr ${#seqfolder[@]} - 1`
		
		for i in `seq 0 $totalfiles`
		do
			
			counts=`expr $counts + 1`
			
			if [ `echo ${files[$i]} | grep -v phix | grep -v "#0.cram" | grep "$filetype$"` ]; then
				fileName=${seqfolder[$i]}/${files[$i]}
				if [ "$counts" == "2" ]; then
					getMetaData $fileName
				else
					getMetaData $fileName | sed -e '1d'
				fi
			fi
		done
	else
		echo "Can not find any data for the above combination"
	fi
done < $inputFile > $title.info
sms "Removing imeta temp files"
rm -f *.imeta


