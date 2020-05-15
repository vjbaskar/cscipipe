
if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args
	-t|--title: title
	-s|--sampleFile: inputSampleFile

	inputSampleFile: <star tab> <name> <condn> <type>
	BCOR_shRNA1.star.ReadsPerGene.out.tab	BCOR_shRNA1	BCOR_shRNA	paired-end
	BCOR_shRNA2.star.ReadsPerGene.out.tab	BCOR_shRNA2	BCOR_shRNA	paired-end
	Control_shRNA1.star.ReadsPerGene.out.tab	Control_shRNA1	Control	paired-end
	Control_shRNA2.star.ReadsPerGene.out.tab	Control_shRNA2	Control	paired-end
	Control2_shRNA1.star.ReadsPerGene.out.tab	Control2_shRNA1	Control2	paired-end
	Control2_shRNA2.star.ReadsPerGene.out.tab	Control2_shRNA2	Control2	paired-end
	Ring1_RNF2_shRNA1.star.ReadsPerGene.out.tab	Ring1_RNF2_shRNA1	Ring1_RNF2_shRNA	paired-end
	Ring1_RNF2_shRNA2.star.ReadsPerGene.out.tab	Ring1_RNF2_shRNA2	Ring1_RNF2_shRNA	paired-end
"
	exit 2
fi


#### Once counting is done, create counts.txt file
sms "Creating count table"
i=0
#t=`tempfile -d ./` 
t=`${mconf_tempfile_cmd}`
while read line
do
	l=($line)
	file=${l[0]}
	name=${l[1]}
	condition=${l[2]} 	
	echo "Reading $file >> $name"
	 if [ $i -eq 0 ]; then
		 echo "gene_name" | tee $title.unstranded.tsv $title.yes.frSS.SOLiD.tsv $title.rev.frFS.Illum.tsv > /dev/null
		 awk ' { print $1 } ' $file | grep -v N_ | tee -a $title.unstranded.tsv $title.yes.frSS.SOLiD.tsv $title.rev.frFS.Illum.tsv > /dev/null
	fi
		# unstranded
		echo $name > $t
		cat $file | grep -v N_ | awk ' { print $2 } ' >> $t
		paste $title.unstranded.tsv $t > $t.us
		mv $t.us $title.unstranded.tsv

		echo $name > $t
		cat $file | grep -v N_ | awk ' { print $3 } ' >> $t
		paste $title.yes.frSS.SOLiD.tsv $t > $t.ss
		mv $t.ss $title.yes.frSS.SOLiD.tsv
		
		echo $name > $t
		cat $file | grep -v N_ | awk ' { print $4 } ' >> $t
		paste $title.rev.frFS.Illum.tsv $t > $t.fs
		mv $t.fs $title.rev.frFS.Illum.tsv
	
	i=`expr $i + 1`
done < $sampleFile
rm $t
