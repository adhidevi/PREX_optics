#!/bin/bash
rm -f TextFiles/Qsq_edgeCut_allRunsL.csv
rm -f TextFiles/Qsq_edgeCut_allRunsR.csv
while read date target runL runR Ebeam
do
   echo $date $target $runL $runR $Ebeam
   rm -f TextFiles/Qsq_pCut_run$runL.csv
   rm -f TextFiles/Qsq_pCut_run$runR.csv
   root -l -b -q "Qsq_only_pCut.C($runL,$Ebeam)"
   root -l -b -q "Qsq_only_pCut.C($runR,$Ebeam)"
done<TextFiles/good_run.list
pdfunite ./plots/Qsq_pCut_run*.pdf ./plots/all_Qsq_pCut.pdf
rm -f ./plots/Qsq_pCut_*
