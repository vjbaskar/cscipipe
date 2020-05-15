# All args come from cscipipe command

Help() { 
echo -en "\033[31m"
cat << EOM
++ Mandatory args: 
-t|--title
-c|--comments
-f|--inputFile: sgRNA summary file from mageck
-v|--qval: [ 0 - 1], eg. 0.2
-0|--cmd0: total sgRNA/gene > cmd0 with fdr <= qval, eg. 3
-1|--cmd1: list of genes interested.

EOM
echo -en "\033[30m"
}


if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile qval cmd0`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

$mconf_R -e "rmarkdown::render('$mconf_installdir/rmds/mageckXYsgRNA.Rmd',output_dir='$cwd/', knit_root_dir = '$cwd', output_file='$cwd/$title.$comments.xy.html')"  

