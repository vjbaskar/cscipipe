
if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: title, gtf, countsFile
	++ recommended: mem=10000   procs=1
	
	organism options: hs,mm

	***** inputFile format
	bamfile1.bam control1.bam
	bamfile2.bam control1.bam
"
	exit 2
fi


#### deseqFPKM.R constructs count tables and then runs FPKM calculation using DESeq
#echo "Rscript $R_SCRIPTS_DIR/deseqFPKM.R $gtf $inputFile
#Rscript $R_SCRIPTS_DIR/plotDESeq.R" > $inputFile.init.sh
#$SCRIPT_DIR/farm..sub -q small -n fpkm -m 4000 -p $THREADS -s $inputFile -j htseqCount -c "sh $inputFile.init.sh" 
#$SCRIPT_DIR/farm..sub -q small -n fpkm -m 4000 -p $THREADS -s $inputFile -c "sh $inputFile.init.sh" -j htseqCount
countsfile=$1; shift
gtf=$1; shift

echo $countsfile
echo $gtf

#### Run pairwise-deseq
sms "Pairwise DEX analysis"
echo -n "" > $inputFile.pw.commands
while read line
	do
		    l=($line)
			c1=${l[0]}
		    c2=${l[1]}
		    # Rscript $R_SCRIPTS_DIR/deseq.R $gtf counts.txt $inputFile $c1 $c2 > ${c1}_VS_$c2.log 2>&1
		   echo "Rscript $mconf_installdir/deseq1.R $gtf counts.txt $inputFile $c1 $c2 > ${c1}_VS_$c2.log 2>&1" >> $inputFile.pw.commands
		   echo "Rscript $mconf_installdir/deseq2.R $gtf counts.txt $inputFile $c1 $c2 > ${c1}_VS_$c2.log 2>&1" >> $inputFile.pw.commands
	done < $sampleFile

#### Submit comparisons as job array
#$SCRIPT_DIR/farm..sub -q small -s deseq -m $MEM -p $THREADS -n $inputFile -j fpkm.$inputFile -a $inputFile.pw.commands -t 20 
echo "$SCRIPT_DIR/farm..sub -q small -s deseq -m $MEM -p $THREADS -n $title -a $title.pw.commands -t 20"
