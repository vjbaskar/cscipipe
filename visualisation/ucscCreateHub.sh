
if [ $args_return  -gt 0 ]; then
print_stderr "
	++ Mandatory args:
	
	options: 
	-t|--title = short name of the hub
	-c|--comments = long name of the hub
	-o|--organism = ucsc organism names such as hg19,hg38,mm10 etc
	-O|--outpt = output directory in web-bfint (your webserver) eg. /data/scratch/project/vm11/GONIA/2019.02.Gonia.HFMX_DNMT3AR882H_FLT3.RRoxBS
	-f|--inputFile
	
	inputFile format:
	-----------------
	file1.bw sample1	cond1	10,12,255
	file2.bw sample2	cond1	255,40,1
	
	Get rgb colours from https://www.rapidtables.com/web/color/RGB_Color.html
"
	exit 2
fi

if [ $args_return  -gt 0 ]; then
	Help
	exit 2
fi

mandatory_fails=`mandatory title comments inputFile organism outpt`
if [ `echo "$mandatory_fails" | wc -w` -gt 0 ]; then
	warnsms "mandatory check failed"
	errorsms "Set: $mandatory_fails"
fi

p=`pwd`
BASE_DIR=$mconf_installdir
name=$title
genome=$organism
webaddr="web-bfint-01"
webscratch=$outpt
inputFile=`readlink -f $inputFile`
#### generate bigwig track
getbw() {
	file=$1
        sample=$2
        cond=$3
        colour=$4
	bw=$file
	track_name=`echo $bw | sed -e 's/.bw$//g' | sed -e 's/.sormadup//g' | sed -e 's/.sormadup//g' `
	c=$colour
	echo "
track $track_name
bigDataUrl $bw
shortLabel $track_name
longLabel $track_name
maxHeightPixel 40:40:11 
color $c
visibility full
smoothingWindow 2
windowingFunction mean
type bigWig 0 1000
autoScale on
viewLimits 0:20
parent $cond
	"
#echo "track type=bigWig name=\"$track_name\" description="" bigDataUrl=$htmlLink/$genome/$file" > /dev/stderr

}


getbb() {
	bb=$1
	c=$2
	track_name=`basename $bb .bb`
	#c=`colour`
	echo "
track $track_name
bigDataUrl $bb
shortLabel $track_name
longLabel $track_name
maxHeightPixel 40:40:11 
color $c
type bigBed
visibility dense
viewLimits 0:20
	"

echo "track type=bigNarrowPeak name=\"$track_name\" description="" bigDataUrl=$htmlLink/$genome/$file" > /dev/stderr
echo "track type=bigBed name=\"$track_name\" description="" bigDataUrl=$htmlLink/$genome/$file" > /dev/stderr
}

createSuperTrack() {
	name=$1
	echo "
track $name
shortLabel $name
longLabel $name
superTrack on show
"

}

createAggTrack() {
	sample=$1
	cond=$2
echo "
track $sample
type bigWig
container multiWig
shortLabel $sample
longLabel $sample
visibility full
aggregate transparentOverlay
showSubtrackColorOnUi on
maxHeightPixels 500:100:8
parent $cond on
"

}

#### generate hubs
#longLabel=`pwd  | xargs basename`
#cat > hub.txt <<'EOF' 
echo "hub $title
shortLabel $title
longLabel $comments
genomesFile genomes.txt
email vm11@sanger.ac.uk
descriptionUrl desc.html" > hub.txt

#### trackDb file

cat <<EOF > genomes.txt
genome $genome
trackDb $genome/trackDb.txt
EOF

#### Copy files to genome

mkdir -p $genome
cp *.bw $genome

#### Populate the trackDb

htmlLink=`echo $webscratch | sed -e 's;/data;ftp://ngs.sanger.ac.uk/;g'`

cd $genome
echo -n "" > trackDb.txt

conds=`cat $inputFile | awk ' { print $3 } ' | sort -u`
for c in $conds
do
	createSuperTrack $c
done >> trackDb.txt


tmpfile=`$mconf_tempfile_cmd`
cat $inputFile | awk ' { print $2,$3 } ' | sort -u > $tmpfile
while read c
do
	temp=($c)
	sample=${temp[0]}
	cond=${temp[1]}
#	createAggTrack $sample $cond
done < $tmpfile >> trackDb.txt

while read line
do
	l=($line)
	file=${l[0]}
	sample=${l[1]}
	cond=${l[2]}
	colour=${l[3]}
	
	suffix=`echo $file | awk -F "." ' { print $NF } '`
	if [ "$suffix" == "bw" ]; then
		sms "Processing as bigWig >> $file"
		getbw $file $sample $cond $colour 
		#cp $cwd/$file .
	fi		
	if [ "$suffix" == "bb" ]; then
		sms "Processing as bigBed >> $file"
		getbb $file $colour $sample 
		cp $cwd/$file .
	fi		
	
	
done < $inputFile >> trackDb.txt
rm -f $tmpfile

cd $cwd
scp -rp $organism hub.txt genomes.txt $webaddr:$webscratch/
echo "Your html hub file is is ${htmlLink}/hub.txt"
