#!/bin/bash
while read date target run
do
   echo $date $target $run
   rm -f TextFiles/output_pCut_run$run.csv
   root -l -b -q "AsymForAccPcut.C($run,\"$date\",\"$target\")"
done<TextFiles/good_runR.list
pdfunite ./plots/all_pCutR_run*.pdf ./plots/all_pCut_runs_RHRS.pdf
