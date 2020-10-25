#!/bin/bash
while read date target run
do
   echo $date $target $run
   rm -f TextFiles/output_pCut_run$run.csv
   while read cut
   do
	echo $cut
#	root -l -b -q "AsymForAccXcut.C($run,$cut)"
	root -l -b -q "AsymForAccPcut.C($run,$cut)"
#	root -l -b -q "AsymForAccYcutLo.C($run,$cut)"
#	root -l -b -q "AsymForAccYcutHi.C($run,$cut)"
   done<TextFiles/pCut.list
#   done<TextFiles/pCutL.list
   pdfunite ./plots/all_pCut_*_run$run.pdf ./plots/plots_pCut_run$run.pdf
   rm -f ./plots/all*
done<TextFiles/good_runR.list
