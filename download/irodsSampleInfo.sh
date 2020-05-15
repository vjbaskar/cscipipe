# Gets sample information given a set of sample ids
# All args come from cscipipe command



Help() { 
echo -en "\033[31m"
cat << EOM

++ Mandatory args: 
	-t|--title
	-f|--inputFile	

	***** inputFile format
	4010STDY7036858
	4010STDY7036859
	4010STDY7036860
	4010STDY7036861
	4010STDY7036862
	4010STDY7036863
	4010STDY7036864


++ Packages used
irods tools


EOM
echo -en "\033[30m"
}



if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title inputFile`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi




function getMetaData(){
#getMetaData /seq/21908/21908_1#1.cram
    fileName=$1
    t=`tempfile -d ./`
    imeta ls -d $fileName > $t
    #errorsms `cat $t`

    fn=`cat $t | grep -A 1 "sample_supplier_name" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    lane=`cat $t | grep -A 1 "lane" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    tagIndex=`cat $t | grep -A 1 "tag_index" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    sample_supplier_name=`cat $t | grep -A 1 "sample_supplier_name" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    seqchksum=`cat $t | grep -A 1 "seqchksum" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    total_reads=`cat $t | grep -A 1 "total_reads" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    library_type=`cat $t | grep -A 1 "library_type" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    sample_donor_id=`cat $t | grep -A 1 "sample_donor_id" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    reference=`cat $t | grep -A 1 "reference" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    id_run=`cat $t | grep -A 1 "id_run" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    is_paired_read=`cat $t | grep -A 1 "is_paired_read" | awk -F ":" ' $1 ~ "value"  { print $2 } ' | sed -e 's/ //g'`
    printf "sample_id\tfile_name\tid_run\tsample_donor_id\tsample_supplier_name\treference\tlibrary_type\tis_paired_read\n"
    printf "${sample}\t$fileName\t$id_run\t$sample_donor_id\t$sample_supplier_name\t$reference\t$library_type\t$is_paired_read\n"
    rm -f $t
}

filetype="cram"
rm -f $title.info
counts=0
while read line
do
	
	l=($line)
	sample=${l[0]}

	sms "Retrieving data from irods <=> $sample " 

	ofile=${sample}
	targetFolder=$sample
	rm -rf targetFolder
	imetafile=${ofile}.imeta
	#imeta qu -z seq -d id_run = $id_run and target = 1 and manual_qc = 1 and sample = $sample > $imetafile
	#imeta qu -z seq -d  target = 1 and manual_qc = 1 and sample = $sample > $imetafile
	imeta qu -z seq -d  sample = $sample > $imetafile
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
			if [ `echo ${files[$i]} | grep -v phix | grep -v "#0.cram" | grep "$filetype$"` ]; then
				fileName=${seqfolder[$i]}/${files[$i]}
				counts=`expr $counts + 1`
				if [ "$counts" == "1" ]; then	
					getMetaData $fileName >> $title.info
					
				else
					getMetaData $fileName | sed -e '1d' >> $title.info
				fi
			fi
		done
	else
		errorsms "Can not find any data for the above combination"
	fi
done < $inputFile 
sms "Removing imeta temp files"
rm -f *.imeta



