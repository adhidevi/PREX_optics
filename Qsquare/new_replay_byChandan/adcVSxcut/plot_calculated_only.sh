#!/bin/bash
rm -f ./TextFiles/adc_x_p_cut_cpmpare.csv
rm -f ./temp1/*.pdf
IFS=$'\t'
while read date target run Ebeam
do
	echo $date $target $run $Ebeam
        root -l -b -q "Qsq_only_x_p_adc_cut.C($run,\"$target\",\"$date\",$Ebeam)"
done < ./TextFiles/good_runs.list
pdfunite ./temp1/*.pdf ./plots/LR_qsq_x_p_adc_cut.pdf
rm -rf ./temp1/*
