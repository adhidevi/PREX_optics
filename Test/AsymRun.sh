#!/bin/bash
while read date target run 
do
	echo $date $target $run
	rm -f ./TextFiles/output_adcCut_run$run.csv
	rm -f ./plots/plots_adcCut_run$run.pdf
	root -l -b -q "AsymForAccCut.C($run,\"$date\",\"$target\")"
done<TextFiles/q2runsR.list
pdfunite ./temp/plots_* ./plots/Asym_plots_allRuns.pdf
rm -f ./temp/Asym_plots*
