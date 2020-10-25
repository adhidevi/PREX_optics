#!/bin/sh
rm -f out_fileL.txt
rm -f out_fileR.txt
rm -f ./temp/*.pdf
run_list="run_list3.list"
IFS=$'\t'
while read date target runL runR
do
	echo "plotting for: " $date $target $runL
	root -l -b -q "plotQsqOnlyUS.C($runL,\"$target\",\"$date\")"
	echo "plotting for: " $date $target $runR
	root -l -b -q "plotQsqOnlyUS.C($runR,\"$target\",\"$date\")"
done < $run_list
pdfunite ./temp/*.pdf ./plots/Qsquare_LR.pdf
