# Summarise samples
if [ "$step" == "summarise" ]; then
    echo "
    summarize_miso --summarize-samples ${output_dir}/${control_name}/ ${output_dir}/${control_name}/
    summarize_miso --summarize-samples ${output_dir}/${test_name}/ ${output_dir}/${test_name}/
    " > ${output_dir}/summarise.sh
    farm..sub -n summarise -s $runName -c "sh ${output_dir}/summarise.sh" -w -q normal
fi

# Compare samples
if [ "$step" == "compare" ]; then
    echo "compare_miso --compare-samples ${output_dir}/${control_name}/ ${output_dir}/${test_name}/ ${output_dir}/" > ${output_dir}/compare.sh
    farm..sub -n compare -s $runName -c "sh ${output_dir}/compare.sh" -w -q normal
fi

# filter events
if [ "$step" == "filter" ]; then
    #miso_bf="/lustre/scratch117/casm/team163/vm11/KOSTAS/2017.07.NPM1c-Mettl16_RNASeq/MISO/A3SS/AF9_comparisons/AF9_WT_vs_AF9_MET16/bayes-factors/AF9_WT_vs_AF9_MET16.miso_bf"
    miso_bf="${output_dir}/${control_name}_vs_${test_name}/bayes-factors/${control_name}_vs_${test_name}.miso_bf"
    echo "filter_events --filter $miso_bf --num-inc 1 --num-exc 1 --num-sum-inc-exc 10 --delta-psi 0.10 --bayes-factor 10 --output-dir ${output_dir}/filtered/" > ${output_dir}/filter.sh
    farm..sub -n filter -s $runName -c "sh ${output_dir}/filter.sh" -q small
    tar -zcvf ${output_dir}/${control_name}.tar.gz ${output_dir}/${control_name}/
    tar -zcvf ${output_dir}/${test_name}.tar.gz ${output_dir}/${test_name}/
fi

if [ "$step" == "annotate" ]; then
# Rscript miso_analysis.R /Users/vm11/nfs_vijay/GENOMES/MISO/GRCm38/all_events/ NPM1_miso/filtered/NPM1_WT_vs_NPM1_MET16.miso_bf.filtered
    Rscript $R_SCRIPTS/miso_analysis.R ${indexed_events} ${output_dir}/filtered/${control_name}_vs_${test_name}.miso_bf.filtered
fi

if [ "$step" == "log" ]; then
    RED='\033[0;31m'
    NC='\033[0m' # No Color
    total_pass=0
    total_fail=0
    mkdir -p mylogs/success mylogs/fail
    for file in `ls ${output_dir}/*/scripts_output/*`
    do
        success=`grep Success $file | wc -l`
        f=`echo $file | sed -e 's:/:____:g'`
        if [ "$success" == 1 ]; then
            ((total_pass++))
    #       cp -p $file mylogs/success/$f
            echo "PASS: $file --> $f"
        else
            ((total_fail++))
            cp -p $file mylogs/fail/$f
            echo "${RED}FAIL: $file --> $f ${NC}"
        fi
    done
    echo "Log files"
    echo "Total pass = $total_pass"
    echo -e "${RED}Total fail = $total_fail ${NC}"
fi



# farm submission: commandfile var = commandFile
# if [ ! -z "$farm" ]; then
# 	export fileOfCommands=$commandFile
# 	$mconf_installdir/bin/farmsub.sh
# else 
# 	warnsms "Not running as farm: Not recommended"
# 	warnsms "$mconf_bashexec $commandFile"
# fi
