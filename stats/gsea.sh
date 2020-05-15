# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments 
-f|--inputFile = Deseq2 output file.
-0|--cmd0 = cols for getting genenames and stats in deseq. eg 1,5
-1|--cmd1 = number of permutations [ 1000 ]
Outpt: 

++ recommended: 
-F|--farm = flag
-m|--mem=4000   # becomes bigger with bam file size
-p|--procs=1
-q|--queue=small

	
EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile cmd0 cmd1`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi



fields=$cmd0
outputFolder=${title}_gsea
nperm=$cmd1
rankFile=`echo $inputFile | sed -e "s/.txt\|.csv\|.tsv/.gseaInput.rnk/g"`
cut -f $fields $inputFile > $rankFile
echo "Correcting input rank file"
tmpfile=`tempfile -d ./`
tr '[a-z]' '[A-Z]' < $rankFile | grep -v "GENE" | grep -v "BASEMEAN" > $tmpfile
sort -nr -k 2 $tmpfile > $tmpfile.1
echo -e "#gene_name\tstat" > $rankFile
cat $tmpfile.1 >> $rankFile
rm -f $tmpfile*

echo "Running GSEA"

commandFile="${title}.${comments}.cmds"
echo -n "" > $commandFile
for set in  c2.all c4.all c5.all c6.all c7.all c2.cp c2.cp.kegg c2.cgp c5.bp 
do	
	mkdir -p $outputFolder/${set}
	cmd="java -cp $mconf_gseaJarFile -Xmx${mem}m xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/${set}.v5.2.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm $nperm -rnk $rankFile -scoring_scheme classic -rpt_label $comment -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 500 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out $outputFolder/${set} -gui false ; zip ${outputFolder}.${set}.zip $outputFolder/${set}/*"
	echo $cmd >> $commandFile
done

mkdir -p $outputFolder/all
cmd="java -cp $mconf_gseaJarFile -Xmx${mem}m xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.v5.2.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.v5.2.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c4.all.v5.2.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.bp.v5.2.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c7.all.v5.2.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm $nperm -rnk $rankFile -scoring_scheme classic -rpt_label $comment -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 500 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out $outputFolder/all -gui false ; zip $outputFolder/all.zip $outputFolder/all/*"
echo $cmd >> $commandFile

# farm submission: commandfile var = commandFile
if [ ! -z "$farm" ]; then
	export fileOfCommands=$commandFile
	$mconf_installdir/bin/farmsub.sh
else 
	warnsms "Not running as farm: Not recommended"
	$mconf_bashexec $commandFile
fi
