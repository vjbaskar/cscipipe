Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-f|--inputFile
-p|--procs

***** inputFile format [ena text download from their website for a study]
study_accession	sample_accession	secondary_sample_accession	experiment_accession	run_accession	tax_id	scientific_name	instrument_model	library_layout	fastq_ftp	fastq_galaxy	submitted_ftp	submitted_galaxy	sra_ftp	sra_galaxy	cram_index_ftp	cram_index_galaxy
PRJNA473835	SAMN09284766	SRS3356892	SRX4143105	SRR7236860	9606	Homo sapiens	Illumina HiSeq 4000	PAIRED	ftp.sra.ebi.ac.uk/vol1/fastq/SRR723/000/SRR7236860/SRR7236860_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR723/000/SRR7236860/SRR7236860_2.fastq.gz	ftp.sra.ebi.ac.uk/vol1/fastq/SRR723/000/SRR7236860/SRR7236860_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR723/000/SRR7236860/SRR7236860_2.fastq.gz		ftp.sra.ebi.ac.uk/vol1/srr/SRR723/000/SRR7236860	ftp.sra.ebi.ac.uk/vol1/srr/SRR723/000/SRR7236860		
PRJNA473835	SAMN09284766	SRS3356892	SRX4143105	SRR7236861	9606	Homo sapiens	Illumina HiSeq 4000	PAIRED	ftp.sra.ebi.ac.uk/vol1/fastq/SRR723/001/SRR7236861/SRR7236861_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR723/001/SRR7236861/SRR7236861_2.fastq.gz	ftp.sra.ebi.ac.uk/vol1/fastq/SRR723/001/SRR7236861/SRR7236861_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR723/001/SRR7236861/SRR7236861_2.fastq.gz		ftp.sra.ebi.ac.uk/vol1/srr/SRR723/001/SRR7236861	ftp.sra.ebi.ac.uk/vol1/srr/SRR723/001/SRR7236861

**** Packages used
wget

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

downloadData() {
	ftp=$1
	#echo $ftp
	x=`echo "$ftp" | awk -F "\t" '  { print $10 } ' | grep -v "^$" | sed -e 's/;/\n/g'`
	#echo "$x"
	echo "$x" | parallel --eta --jobs $procs wget -nv 
}


while read line
do		
		downloadData "$line"
done < <(cat $inputFile | grep -v "study_accession")
