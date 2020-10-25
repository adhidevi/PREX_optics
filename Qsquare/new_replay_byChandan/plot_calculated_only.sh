#!/bin/bash
rm -f ./TextFiles/theta_LR.csv
rm -f ./temp/*.pdf
IFS=$'\t'
while read date target runL runR Ebeam
do
	echo $date $target $runL $runR $Ebeam
        root -l -b -q "calculated_only_theta.C($runL,$runR,\"$target\",\"$date\",$Ebeam)"
done < ./TextFiles/good_runs.list
pdfunite ./temp/*.pdf ./plots/LR_calculated_Theta.pdf
rm -rf ./temp/*
