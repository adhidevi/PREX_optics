#!/bin/bash
rm -f ./TextFiles/qsquare_diagnostic_quantitiesL.csv
rm -f ./TextFiles/qsquare_diagnostic_quantitiesR.csv
IFS=$'\t'
while read date target run Ebeam
do
	echo $date $target $run $Ebeam
        root -l -b -q "projXY_and_findQsq.C($run,\"$target\",\"$date\",$Ebeam)"
done < ./TextFiles/good_runs.list
pdfunite ./temp1/*.pdf ./plots/qsq_diagnostic_plots.pdf
rm -rf ./temp1/*
