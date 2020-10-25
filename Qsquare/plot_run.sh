#!/bin/bash
rm -f out_fileL.txt
rm -f out_fileR.txt
rm -f ./TextFiles/qsq_analyzerR.csv
rm -f ./TextFiles/qsq_analyzerL.csv
rm -f ./TextFiles/qsq_calculatedR.csv
rm -f ./TextFiles/qsq_calculatedL.csv
rm -f ./temp/*.pdf
IFS=$'\t'
#
while read date target runL runR
do
	echo $date $target $runL $runR
#	root -l -b -q plotThetaPhi.C\($runL,\"$target\"\)
#	root -l -b -q plotThetaPhi.C\($runR,\"$target\"\)
#	root -l -b -q theta.C\($runL,\"$target\"\)
#	root -l -b -q theta.C\($runR,\"$target\"\)
#	root -l -b -q "plotMomentum.C($runL,\"$target\")"
#	root -l -b -q "plotMomentum.C($runR,\"$target\")"
#        root -l -b -q "plotQsqOnlyUS.C($runL,\"$target\",\"$date\")"
#        root -l -b -q "plotQsqOnlyUS.C($runR,\"$target\",\"$date\")"
        root -l -b -q "Qsq_ana_cal.C($runL,\"$target\",\"$date\")"
        root -l -b -q "Qsq_ana_cal.C($runR,\"$target\",\"$date\")"
done < ./TextFiles/good_runs.list
pdfunite ./temp/*.pdf ./plots/important_plots_Qsq_runs.pdf
rm -rf ./temp/*
