#!/bin/sh


if [ $args_return -gt 0 ]; then
print_stderr "
	++ Mandatory args: 
	-t|--title
	-c|--comments
	-f|--inputFile
	-o|--organism : star index
	-F|--farm
	-P|--progArgs : additional prog args (eg: --limitOutSJcollapsed 3200000)
	++ recommended: 
	-m|--mem=32000   
	-p|--procs=1
	-r|--paired: 0/1
	
	++ optional:

	
	inputFile: 
	SRR6144783	WT1
	SRR6144784	WT2
"
	exit 2
fi

#threads=$1; shift
#paired=$1; shift
#genome=$1; shift
#namefile=$1; shift
#index=/nfs/users/nfs_v/vm11/GENOMES/STAR_INDEX/${genome}
#mem=32000


threads=$procs
genome=$organism
namefile=$inputFile
index=$mconf_starindex/${organism}
index=${organism}


seAlign() {
	id=$1
	name=$2
	STAR_INPUT="$id.fq.gz"
	OUTPUT_PREFIX="$name.star."
	echo "
		STAR --runThreadN  ${threads}  ${progArgs} \
		--genomeDir ${index} \
		--readFilesIn ${STAR_INPUT} \
		--outFileNamePrefix ${OUTPUT_PREFIX} \
		--readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
		--outSAMunmapped Within \
		--quantMode GeneCounts
		$mconf_biobambam2/bamsormadup threads=$threads SO=coordinate M=temp.metrics indexfilename=${OUTPUT_PREFIX}sormadup.bam.bai < ${OUTPUT_PREFIX}Aligned.out.bam > ${OUTPUT_PREFIX}sormadup.bam
	"
}


peAlign() {
	id=$1
	name=$2
	STAR_INPUT="$id.pair1.fq.gz $id.pair2.fq.gz"
	OUTPUT_PREFIX="$name.star."
	echo "  set -euxo pipefail
		STAR --runThreadN  ${threads} ${progArgs} \
		--genomeDir ${index} \
		--readFilesIn ${STAR_INPUT} \
		--outFileNamePrefix ${OUTPUT_PREFIX} \
		--readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
		--outSAMunmapped Within \
		--quantMode GeneCounts 
		$mconf_biobambam2/bamsormadup threads=$threads SO=coordinate M=temp.metrics indexfilename=${OUTPUT_PREFIX}sormadup.bam.bai < ${OUTPUT_PREFIX}Aligned.out.bam > ${OUTPUT_PREFIX}sormadup.bam
	"
}

echo -n "" > $title.star.cmds
while read line
do
		l=($line)
		id=${l[0]}
		name=${l[1]}

		if [ $paired -eq 1 ]; then
			echo "Running in paired end mode:	$id <> $name"
			peAlign $id $name > $id.$name.star.sh
		else
			echo "Running in single end mode:	$id <> $name"
			seAlign $id $name > $id.$name.star.sh
		fi
		#sub_ncores normal star $mem $threads $name "sh $id.$name.star.sh"
		echo "sh $id.$name.star.sh" >> $title.star.cmds

done < $namefile


# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export fileOfCommands=$title.star.cmds
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	warnsms "$mconf_bashexec $commandFile"
fi


