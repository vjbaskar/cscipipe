
if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args: organism, paired, qval, inputFile, mem, procs
	++ recommended: mem=10000   procs=1
	
	organism options: hs,mm

	***** inputFile format
	bamfile1.bam control1.bam
	bamfile2.bam control1.bam

"
	exit 2
fi

echo $file1
echo $file0
