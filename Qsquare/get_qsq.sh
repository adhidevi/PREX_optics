#!/bin/sh
rm -f ./temp/*.pdf
run_list="good_runs.list"
IFS=$'\t'
while read date target runL runR
do
	echo "plotting for: " $date $target $runL
	root -l -b -q "Qsq_ana_cal.C($runL,\"$target\")"
	echo "plotting for: " $date $target $runR
	root -l -b -q "Qsq_ana_cal.C($runR,\"$target\")"
done < $run_list
pdfunite ./temp/*.pdf ./plots/Qsquare_compare_ana_cal.pdf
