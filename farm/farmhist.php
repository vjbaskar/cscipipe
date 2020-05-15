#!/usr/bin/env php
<?php
function jobData($filename, $p, $f) {

	$fileContents = file($filename);
	for($i=0;$i<=count($fileContents);$i++){
		$line=$fileContents[$i];
		$line=chop($line,"\n");
		$t = preg_split("/\s+/",$line);
		if($t[0] == "Subject:"){  $id = str_replace(":","",$t[2]);   }
		if(preg_match("/^Job </",$line)){  $jobName= preg_replace("[\<|\>]","",$t[1]); }
		#$p="/".$p."/";
		#echo "$p" ; 
		if(preg_match("/CPU/",$line)){  $cputime= $t[4]; }
		if(preg_match("/Max Memory/",$line)){  $maxmem= $t[4]; }
		if(preg_match("/Average Memory/",$line)){  $avgmem= $t[4]; }
		if(preg_match("/Total Requested Memory/",$line)){  $totalmem= $t[5]; }
		if(preg_match("/Max Processes/",$line)){  $maxProcs= $t[4]; }
		if(preg_match("/Started at/",$line)){  $startTime= preg_replace("/Started at /","",$line); }
		if(preg_match("/Results reported on /",$line)){  $endTime= preg_replace("/Results reported on /","",$line); }
		if(preg_match("/ was used as the working directory./",$line)){  $cwd = preg_replace("/ was used as the working directory./","",$line); }
		if(preg_match("/Resource usage summary/",$line)){  $status= chop($fileContents[$i-2],"\n"); }
		if(preg_match("/The output/",$line)){  break; }

	}
	if ($f == 1) {
		echo "status\tfilename\tid\tjobName\t<startTime>\t<endTime>\tcputime\tmaxmem\tavgmem\ttotalmem\tmaxProcs\tcwd\n";
		echo "------\t--------\t--\t-------\t-----------\t---------\t-------\t------\t------\t--------\t--------\t---\n";
	}
	if (preg_match("/$p/", $jobName)) {
		echo "$status\t$filename\t$id\t$jobName\t<$startTime>\t<$endTime>\t$cputime\t$maxmem\t$avgmem\t$totalmem\t$maxProcs\t$cwd\n";
	} elseif ($p == "all"){
		echo "$status\t$filename\t$id\t$jobName\t<$startTime>\t<$endTime>\t$cputime\t$maxmem\t$avgmem\t$totalmem\t$maxProcs\t$cwd\n";
	}
}

$progArgs = getEnv("progArgs");

if ($progArgs == "") {
	echo "---------\n";
	echo "Required args: ";
	echo "-P|--progArgs [ all | parts of job names, eg: test, bwa ... etc ]" ;
	echo "\n\n";
	exit(2);

}

$dirhandle = opendir(".bsub/");
$f=0;
while (($file = readdir($dirhandle)) !== false){
  if(preg_match("/.farm$/",$file)){
	$file=".bsub/".$file;
	$f = $f + 1;
	jobData($file, $progArgs, $f);
  }
}
closedir($dirhandle);

?>